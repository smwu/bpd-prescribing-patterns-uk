#===============================================================================
# Create Sankey plots and heatmaps of Patterns of treatment persistence and 
# medication class switching within one year of first-recorded BPD diagnosis
# 
# Data:
#   - Data/bipolar_cohort_ads.rdata: Antidepressant prescriptions
#   - Data/bipolar_cohort_moodstabs.rdata: Mood stabiliser prescriptions
#   - Data/bipolar_cohort_aps.rdata: Antipsychotic prescriptions
#   - Data/bipolar_cohort.rdata: Patient information 
# 
# Outputs:
#   - Shiny_App/data/trajectories_limit_counts.rds: Treatment trajectories data for Shiny app 
#   - Outputs/overall_switches_alluvial: Sankey plot of overall treatment patterns within 1 year of BPD diagnosis
#   - Outputs/overall_switches_alluvial_first_state: Sankey plot of overall treatment patterns, colored by first state 
#   - Outputs/overall_switches_alluvial_60day: Sankey plot of overall treatment patterns using 60-day treatment window
#   - Outputs/overall_switches_alluvial_1day: Sankey plot of overall treatment patterns using 1-day treatment window
#   - Shiny_App/data/trajectories_limit_li_counts.rds: Lithium-specific treatment trajectories data for Shiny app 
#   - Outputs/overall_switches_alluvial_lithium: Sankey plot of lithium-focused treatment patterns within 1 year of BPD diagnosis
#   - Outputs/overall_switches_alluvial_lithium_first_state: Sankey plot of lithium-focused treatment patterns, colored by first state 
#   - Outputs/overall_switches_alluvial_AD: Sankey plot of treatment patterns, subsetted to antidepressant initial state
#   - Outputs/initial_transition_heatmap: Heatmap of transition from first- to second-line treatment
#   - Outputs/initial_transition_li_heatmap: Lithium-focused heatmap of transition from first- to second-line treatment
#   - Outputs/second_transition_heatmap: Heatmap of transition from second- to third-line treatment
#   - Outputs/second_transition_li_heatmap: Lithium-focused heatmap of transition from second- to third-line treatment
# ==============================================================================


# Clear memory
rm(list = ls())

library(lubridate)
library(gtsummary)
library(tidyverse)
library(ggforce)
library(purrr)
library(ggalluvial)
library(gt)
library(tidylog)
library(data.table)
library(scales)
library(ggrepel)
library(dplyr)
library(ggplot2)
library(scales)
library(ggtext)

# Load data

load("Data/bipolar_cohort_ads.rdata")
load("Data/bipolar_cohort_moodstabs.rdata")
load("Data/bipolar_cohort_aps.rdata")
load("Data/bipolar_cohort.rdata")

wd <- getwd()

# Parameters
overlap_days <- 30  # drugs within 30 days are considered combination
interrupt_cutoff_intervals <- 3  # >= 3 consecutive intervals
max_switches <- 3 # maximum number of switches for sankey plot

#========== Diagnosis data =====================================================

# Get first bipolar diagnosis date
bipolar_diagnosis_date <- bipolar_cohort %>%
  mutate(bipolar_diag_date = as.Date(first_bipolar_date)) %>%
  select(patid, bipolar_diag_date)

# Censor at earliest of death, lcd, or regendate -- include gender
censor_dates <- bipolar_cohort %>%
  select(patid, gender, deathdate, regenddate, lcd) %>%
  left_join(bipolar_diagnosis_date, by = "patid") %>%
  mutate(
    deathdate = as.Date(deathdate),
    regenddate = as.Date(regenddate),
    lcd = as.Date(lcd),
    bipolar_diag_date = as.Date(bipolar_diag_date),
    earliest_exit = pmin(deathdate, regenddate, lcd, na.rm = TRUE),
    diff_since_diag = earliest_exit - bipolar_diag_date, 
    month_since_diag = ifelse(diff_since_diag >= 0,  # (should all be >= 0)
                              ceiling(diff_since_diag / 30.42),  # after dianosis
                              floor(diff_since_diag / 30.42)),   # before diagnosis
    within_12_months = !is.na(month_since_diag) & abs(month_since_diag) <= 12) %>%
  select(patid, gender, earliest_exit, diff_since_diag, month_since_diag, within_12_months)


#========== Prescription data ==================================================

# Combine valproate variants into "Valproate"
bipolar_cohort_moodstabs <- bipolar_cohort_moodstabs %>%
  mutate(ingredient = case_when(ingredient %in% c("Valproate", "Valproate/Valproic", "Valproic") ~ "Valproate",
                                TRUE ~ ingredient))

# Combine prescriptions from all drug classes
ads <- bipolar_cohort_ads %>%
  select(patid, issuedate, ingredient) %>%
  mutate(drug_type = "AD")

aps <- bipolar_cohort_aps %>%
  select(patid, issuedate, ingredient) %>%
  mutate(drug_type = "AP")

moodstabs <- bipolar_cohort_moodstabs %>%
  select(patid, issuedate, ingredient) %>%
  mutate(drug_type = "MS")

# Combine all prescriptions & remove duplicates
all_prescriptions <- bind_rows(ads, aps, moodstabs) %>%
  distinct()

#========= Combine prescription and diagnosis data =============================

# Join with diagnosis data
prescriptions_with_diag <- all_prescriptions %>%
  left_join(bipolar_diagnosis_date, by = "patid")

# Filter to those with prescription data within one year after diagnosis
prescriptions_filtered_oneyr_after <- prescriptions_with_diag %>%
  mutate(diff_days = as.numeric(issuedate - bipolar_diag_date),
         cutoff_date = bipolar_diag_date + 365) %>%
  filter(diff_days >= 0 & diff_days <= 365)

# Num patients with prescription data within one year after diagnosis: 37,208
length(unique(prescriptions_filtered_oneyr_after$patid))

# Prepare prescriptions data: n = 37208
prescriptions_prep <- prescriptions_filtered_oneyr_after %>%
  mutate(
    issuedate = as.Date(issuedate),
    bipolar_diag_date = as.Date(bipolar_diag_date)) %>%
  # Add in censoring and gender information
  left_join(censor_dates %>% 
              select(patid, gender, earliest_exit), 
            by = "patid")

#=============== Get drug combinations in treatment intervals ==================

# Get treatments intervals
# Example for 30-day intervals:
# Interval 1: first_rx_dat (1st day) to 30th day
# Interval 2: Day 31-60 
# Interval 3: Day 61-90 
# Inputs:
#   prescriptions_prep: df of prepared prescriptions data with columns for 
#     patid, issuedate, ingredient, drug_type, bipolar_diag_date, diff_days, 
#     cutoff_date, gender, earliest_exit
#   overlap_days: Number specifying overlap for calculating combinations
#   interrupt_cutoff_intervals: Number of consecutive cutoff intervals considered a 
#     treatment interruption/discontinuation
# Returns: 
#   prescriptions_combs_final: df of treatment intervals with columns for 
#     patid, int_idx, drug_combination, int_start, int_end, first_rx_date, 
#     earliest_exit, gender, max_int_idx, is_none, run_id, run_length
get_treatment_intervals <- function(prescriptions_prep, overlap_days, 
                                    interrupt_cutoff_intervals) {
  prescriptions_intervals <- prescriptions_prep %>%
    group_by(patid) %>%
    mutate(first_rx_date = min(issuedate, na.rm = TRUE)) %>%
    ungroup() %>%
    # Compute interval index from first Rx (1,2,3,...)
    # Get interval start and end points for each patient and interval
    mutate(int_idx = ceiling(as.numeric(issuedate - first_rx_date + 1) / overlap_days),
           int_start = first_rx_date + overlap_days * (int_idx - 1),
           int_end = int_start + overlap_days - 1) %>%
    # Filter to complete intervals (i.e., int_end comes before the 365 cutoff date)
    filter(int_end < cutoff_date)
  
  prescriptions_combs <- prescriptions_intervals %>%
    # Combine multiple drug types within the same interval
    group_by(patid, int_idx) %>%
    summarise(
      drug_combination = {
        dx <- sort(unique(na.omit(drug_type)))
        if (length(dx) == 0) {
          NA_character_ 
        } else {
          paste(dx, collapse = "; ")
        }
      },
      int_start = first(int_start),
      int_end = first(int_end),
      first_rx_date = first(first_rx_date),
      earliest_exit = first(earliest_exit),
      gender = first(gender)) %>%
    ungroup()
  
  
  # If censored day lies within a time interval, all intervals after are censored. 
  # Filter out those with censored date before int_start
  prescriptions_combs_censor <- prescriptions_combs %>%
    mutate(drug_combination = ifelse(earliest_exit < int_start, "Censored", 
                                     drug_combination)) %>%
    filter(earliest_exit > int_start) %>%
    group_by(patid) %>%
    # Get maximum interval for each patient
    mutate(max_int_idx = max(int_idx, na.rm = TRUE)) %>%
    ungroup()
  
  # Fill in intervals with missing prescription
  prescriptions_combs_fill <- prescriptions_combs_censor %>%
    group_by(patid) %>%
    complete(int_idx = full_seq(1:unique(max_int_idx), 1)) %>%
    mutate(drug_combination = replace_na(drug_combination, "None")) %>%
    ungroup()
  
  # Identify interruption runs (>=3 consecutive "None")
  prescriptions_combs_final <- prescriptions_combs_fill %>% 
    arrange(patid, int_idx) %>%
    group_by(patid) %>%
    mutate(
      is_none = drug_combination == "None",
      # run ID to each consecutive block of identical values of is_none
      run_id = data.table::rleid(is_none), 
      # number of rows in each consecutive block
      run_length = ave(is_none, run_id, FUN = length),
      # Flag interruptions/discontinuations
      drug_combination = case_when(
        is_none & (run_length >= interrupt_cutoff_intervals) ~ "Int/Disc",
        TRUE ~ drug_combination
      )
    )
  
  return(prescriptions_combs_final)
}

# n = 37067 after removing incomplete intervals and censoring
prescriptions_combs_final <- get_treatment_intervals(
  prescriptions_prep = prescriptions_prep, overlap_days = overlap_days, 
  interrupt_cutoff_intervals = interrupt_cutoff_intervals)

#================ Derive treatment trajectories ================================

# Derive treatment trajectories
# Input:
#   prescriptions_combs_final: df output from 'get_treatment_intervals()'
#   max_switches: Number specifying maximum number of switches to consider
#   lithium: Boolean specifying if the trajectories are lithium-focused. If so,
#     column "nonli_combo_raw" will be used to flag switches between non-lithium 
#     combos. This differentiates it from "no switch" for the lithium trajectories.
# Output:
#   trajectories_limit: df with trajectories for each patient, containing 
#     columns for patid, states, num_switches, state_1, state_2,..., state_max_switches
derive_trajectories <- function(prescriptions_combs_final, max_switches,
                                lithium = FALSE) {
  # Function to remove consecutive duplicates in a sequence
  collapse_consecutive <- function(x) x[c(TRUE, x[-1] != x[-length(x)])]
  
  # Keep only non-empty intervals, tidy, order
  df <- prescriptions_combs_final %>%
    filter(drug_combination != "None") %>% # Ignore empty intervals
    mutate(drug_combination = str_squish(drug_combination)) %>% # Remove whitespace
    arrange(patid, int_idx)
  
  # Lithium focus: build per-interval states with injected "Switch Non-Lithium" 
  # to mark when there is a switch from one non-lithium combo to another 
  # non-lithium combo
  if (lithium) {
    # Detect within-non-lithium switching using nonli_combo_raw
    df <- df %>%
      group_by(patid) %>%
      mutate(
        # High-level class at this interval
        hl_now  = drug_combination,
        # Is current/previous interval non-lithium?
        is_nonli_now  = hl_now  == "Non-Lithium",
        is_nonli_prev = dplyr::lag(is_nonli_now, default = FALSE),
        # Did the non-lithium subtype change?
        nonli_change = is_nonli_now & is_nonli_prev &
          !is.na(nonli_combo_raw) & !is.na(dplyr::lag(nonli_combo_raw)) &
          nonli_combo_raw != dplyr::lag(nonli_combo_raw),
        # Cumulative count of within–non-lithium changes
        nonli_switch_id = cumsum(dplyr::coalesce(nonli_change, FALSE)),
        # Internal state that distinguishes each within–non-li change
        state_internal = dplyr::case_when(
          is_nonli_now & nonli_switch_id > 0 ~ paste0("Switch Non-Lithium#", 
                                                      nonli_switch_id),
          TRUE ~ hl_now
        )
      ) %>%
      ungroup()
    
    # List state switches
    trajectories <- df %>%
      group_by(patid) %>%
      summarise(
        states = list(collapse_consecutive(state_internal)),
        num_switches = length(unlist(states)),
        .groups = "drop") 
    
  } else {  # No lithium focus, use original paths
    # List state switches
    trajectories <- df %>%
      group_by(patid) %>%
      summarise(
        states = list(collapse_consecutive(drug_combination)),
        num_switches = length(unlist(states)),
        .groups = "drop") 
  }
  
  
  # Pad or truncate to exactly max_switches number of states
  trajectories_limit <- trajectories %>%
    mutate(
      states_limit = lapply(states, function(s) {
        s <- unique(s)
        if (length(s) >= max_switches) {
          s[1:max_switches]
        } else {
          c(s, rep(NA_character_, max_switches - length(s)))
        }
      })) %>%
    tidyr::unnest_wider(states_limit, names_sep = "_") %>%
    dplyr::rename_with(~ paste0("state_", seq_along(.)), starts_with("states_"))
  
  # Strip internal "#k" tags so downstream code sees clean labels
  clean_hash <- function(x) if (is.null(x)) x else gsub("#\\d+$", "", x)
  state_cols <- grep("^state_\\d+$", names(trajectories_limit), value = TRUE)
  trajectories_limit[state_cols] <- lapply(trajectories_limit[state_cols], clean_hash)
  
  return(trajectories_limit)
}

# n = 37067
trajectories_limit <- derive_trajectories(
  prescriptions_combs_final = prescriptions_combs_final, 
  max_switches = max_switches)


# Aggregate data function for Shiny app
# Minimum number of counts is 5
make_safe_counts <- function(df, k = 5, round_to = 5) {
  df %>%
    mutate(
      state_2 = ifelse(is.na(state_2), "No Switch", state_2),
      state_3 = ifelse(is.na(state_3), "No Switch", state_3)
    ) %>%
    count(state_1, state_2, state_3, name = "n") %>%
    mutate(
      n = ifelse(n < k, 0, round(n / round_to) * round_to)
    ) %>%
    filter(n > 0)
}

# # Save for Shiny app
# n = 36910 after suppressing trajectories with <5 counts
trajectories_limit_counts <- make_safe_counts(trajectories_limit)
# saveRDS(trajectories_limit_counts,
#         paste0("bpd-prescribing-uk-shiny/data/",
#                "trajectories_limit_counts.rds"))


#====================== Examine censoring ======================================

# Since, 365/30 = 12.16, at most 12 intervals per person. This would happen if 
# first prescription was on the day of diagnosis and prescriptions continued 
# for a year after

# Distribution of number of 30-day intervals for all patients, with the first 
# interval starting on the date of first prescription after diagnosis
# smaller numbers will happen if an individual: 
# 1) dies early, 
# 2) starts treatment a long time after diagnosis, or 
# 3) discontinues treatment early on
max_intervals <- prescriptions_combs_final %>% 
  group_by(patid) %>%
  summarise(max_int_idx = first(max_int_idx))
hist(max_intervals$max_int_idx)
table(max_intervals$max_int_idx)


## Time from diagnosis to first prescription
time_to_first_prescript <- prescriptions_prep %>%
  group_by(patid) %>%
  summarise(first_rx_date = min(issuedate, na.rm = TRUE),
            bipolar_diag_date = first(bipolar_diag_date)) %>%
  mutate(time_first_rx = as.numeric(first_rx_date - bipolar_diag_date))

hist(time_to_first_prescript$time_first_rx)
summary(time_to_first_prescript$time_first_rx)


#=============== Create Sankey plot ============================================


create_sankey <- function(trajectories_limit, state_levels = NULL, 
                          state_full_labels = NULL,
                          state_colors = NULL, alpha_input = 0.8, 
                          text_size = 3) {
  
  ## Define factor order and colors
  # Order of states in each step
  if (is.null(state_levels)) {
    state_levels <- c(
      "AD", "AP", "MS",
      "AD; AP", "AD; MS", "AP; MS",
      "AD; AP; MS",
      "Int/Disc",           # abbreviated for Interruption/Discontinuation
      "No Switch"
    )
  }
  
  if (is.null(state_colors)) {
    # Define colors
    state_colors <- c(
      "AD" = "#E41A1C",          # bright red
      "AP" = "#377EB8",          # blue
      "MS" = "#FFD92F",          # yellow
      "AD; AP" = "#984EA3",      # purple
      "AD; MS" = "#FF7F00",      # orange
      "AP; MS" = "#4DAF4A",      # green
      "AD; AP; MS" = "#994F00",  # brown
      "Int/Disc" = "#999999",    # gray
      "No Switch" = "#F5F5F5"
    )
  }
  
  # Map abbreviation states
  if (is.null(state_full_labels)) {
    state_full_labels <- c(
      "AD" = "Antidepressants", 
      "AP" = "Antipsychotics", 
      "MS" = "Mood Stabilisers",
      "AD; AP" = "Antidepressants & Antipsychotics", 
      "AD; MS" = "Antidepressants & Mood Stabilisers", 
      "AP; MS" = "Antipsychotics & Mood Stabilisers",
      "AD; AP; MS" = "Antidepressants & Antipsychotics & Mood Stabilisers",
      "Int/Disc" = "Interruption/Discontinuation",    
      "No Switch" = "No Switch")
  }
  
  
  flow_table_ggalluvial <- trajectories_limit %>%
    filter(!is.na(state_1)) %>%
    mutate(
      state_2 = ifelse(is.na(state_2), "No Switch", state_2),
      state_3 = ifelse(is.na(state_3), "No Switch", state_3)
    )
  # Number in each unique trajectory: use existing counts if found
  if (!"n" %in% names(flow_table_ggalluvial)) {
    flow_table_ggalluvial <- count(flow_table_ggalluvial, 
                                   state_1, state_2, state_3, name = "n")
  } else {
    flow_table_ggalluvial <- flow_table_ggalluvial %>% 
      select(state_1, state_2, state_3, n)
  }
  # Create a unique trajectory ID before pivoting
  flow_table_ggalluvial <- flow_table_ggalluvial %>% 
    mutate(traj_id = row_number()) %>%
    pivot_longer(
      cols = starts_with("state_"),
      names_to = "step",
      values_to = "state"
    )
  # Convert state to factor
  flow_table_ggalluvial <- flow_table_ggalluvial %>%
    mutate(state = factor(state, levels = state_levels))
  
  
  ### Update with percentages
  
  # 1) Make y a per-column proportion and build a per-stratum label
  flow_table_ggalluvial_lab <- flow_table_ggalluvial %>%
    group_by(step) %>%
    mutate(step_total = sum(n),
           pct = n / step_total) %>%                 # used for y
    group_by(step, state) %>%
    mutate(
      stratum_pct = sum(n) / first(step_total),
      label = sprintf("%s (%s)", state, percent(stratum_pct, accuracy = 1))
    ) %>%
    ungroup()
  
  # 2) Plot with y = pct, and show the precomputed label inside each box
  alluvial_plot <- ggplot(
    flow_table_ggalluvial_lab,
    aes(
      x = step,
      stratum = state,
      alluvium = traj_id,
      y = pct,
      fill = state
    )
  ) +
    geom_flow(stat = "alluvium", 
              lode.guidance = "frontback", 
              alpha = alpha_input,
              knot.pos = 0.5) +
    # Wider boxes here (was 0.25)
    geom_stratum(width = 0.40, color = "gray40") +
    geom_text(
      stat = "stratum",
      aes(label = label),
      size = text_size,
      fontface = "bold",
      check_overlap = TRUE,
      color = "black"
    ) +
    scale_x_discrete(
      labels = c("First", "Second", "Third"),
      # Less white space on the left/right
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_fill_manual(values = state_colors,
                      breaks = names(state_full_labels),
                      labels = state_full_labels,
                      name = "Drug Combination") +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      y = "Percent of Patients",
      x = "Treatment Step",
      fill = "Drug Combination"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 13),
      axis.text.x = element_text(size = 11),
      # Slightly tighter outer margins
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
  
  
  alluvial_plot
}


alluvial_plot <- create_sankey(trajectories_limit = trajectories_limit_counts)
alluvial_plot

# # Save the plot
# ggsave(
#   plot = alluvial_plot,
#   filename = paste0(wd, "Outputs/", "overall_switches_alluvial", today(), ".png"),
#   dpi = 300, width = 13, height = 9, bg = "white"
# )

# # Save as PDF
# ggsave(
#   plot = alluvial_plot,
#   filename = paste0(wd, "Outputs/", "overall_switches_alluvial", today(), ".pdf"),
#   dpi = 300, width = 13, height = 9, bg = "white"
# )


# Create Sankey plot, coloring by first state

create_sankey_first_state <- function(trajectories_limit, state_levels = NULL, 
                                      state_full_labels = NULL,
                                      state_colors = NULL, alpha_input = 0.8,
                                      text_size = 3) {
  
  ## Define factor order and colors
  # Order of states in each step
  if (is.null(state_levels)) {
    state_levels <- c(
      "AD", "AP", "MS",
      "AD; AP", "AD; MS", "AP; MS",
      "AD; AP; MS",
      "Int/Disc",           # abbreviated for Interruption/Discontinuation
      "No Switch"
    )
  }
  
  # Define colors
  if (is.null(state_colors)) {
    state_colors <- c(
      "AD" = "#E41A1C",          # bright red
      "AP" = "#377EB8",          # blue
      "MS" = "#FFD92F",          # yellow
      "AD; AP" = "#984EA3",      # purple
      "AD; MS" = "#FF7F00",      # orange
      "AP; MS" = "#4DAF4A",      # green
      "AD; AP; MS" = "#994F00",  # brown
      "Int/Disc" = "#999999",    # gray
      "No Switch" = "#F5F5F5"
    )
  }
  
  # Map abbreviation states
  if (is.null(state_full_labels)) {
    state_full_labels <- c(
      "AD" = "Antidepressants", 
      "AP" = "Antipsychotics", 
      "MS" = "Mood Stabilisers",
      "AD; AP" = "Antidepressants & Antipsychotics", 
      "AD; MS" = "Antidepressants & Mood Stabilisers", 
      "AP; MS" = "Antipsychotics & Mood Stabilisers",
      "AD; AP; MS" = "Antidepressants & Antipsychotics & Mood Stabilisers",
      "Int/Disc" = "Interruption/Discontinuation",    
      "No Switch" = "No Switch")
  }
  
  
  flow_table_ggalluvial <- trajectories_limit %>%
    filter(!is.na(state_1)) %>%
    mutate(
      state_2 = ifelse(is.na(state_2), "No Switch", state_2),
      state_3 = ifelse(is.na(state_3), "No Switch", state_3)
    ) 
  # Number in each unique trajectory: use existing counts if found
  if (!"n" %in% names(flow_table_ggalluvial)) {
    flow_table_ggalluvial <- count(flow_table_ggalluvial, 
                                   state_1, state_2, state_3, name = "n")
  } else {
    flow_table_ggalluvial <- flow_table_ggalluvial %>% 
      select(state_1, state_2, state_3, n)
  }
  # Create a unique trajectory ID before pivoting
  flow_table_ggalluvial <- flow_table_ggalluvial %>% 
    # Create a unique trajectory ID before pivoting
    mutate(traj_id = row_number(),
           # Drug Combination in initial state
           first_state = state_1) %>%
    pivot_longer(
      cols = starts_with("state_"),
      names_to = "step",
      values_to = "state"
    )
  
  flow_table_ggalluvial <- flow_table_ggalluvial %>%
    mutate(state = factor(state, levels = state_levels),
           first_state = factor(first_state, levels = state_levels))
  
  
  ### Update with percentages
  
  # 1) Make y a per-column proportion and build a per-stratum label
  flow_table_ggalluvial_lab <- flow_table_ggalluvial %>%
    group_by(step) %>%
    mutate(step_total = sum(n),
           pct = n / step_total) %>%                 # used for y
    group_by(step, state) %>%
    mutate(
      stratum_pct = sum(n) / first(step_total),
      label = sprintf("%s (%s)", state, percent(stratum_pct, accuracy = 1))
    ) %>%
    ungroup()
  
  # 2) Plot with y = pct, and show the precomputed label inside each box
  alluvial_plot <- ggplot(
    flow_table_ggalluvial_lab,
    aes(
      x = step,
      stratum = state,
      alluvium = traj_id,
      y = pct,
      fill = first_state
    )
  ) +
    geom_flow(stat = "alluvium", 
              lode.guidance = "frontback", 
              alpha = alpha_input,
              knot.pos = 0.5) +
    # Wider boxes here (was 0.25)
    geom_stratum(fill = "gray90", width = 0.40, color = "gray40") +
    geom_text(
      stat = "stratum",
      aes(label = label),
      size = text_size,
      fontface = "bold",
      check_overlap = TRUE,
      color = "black"
    ) +
    scale_x_discrete(
      labels = c("First", "Second", "Third"),
      # Less white space on the left/right
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_fill_manual(values = state_colors,
                      breaks = names(state_full_labels),
                      labels = state_full_labels,
                      name = "Drug Combination") +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      y = "Percent of Patients",
      x = "Treatment Step",
      fill = "Drug Combination"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 13),
      axis.text.x = element_text(size = 11),
      # Slightly tighter outer margins
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
  
  alluvial_plot
}


alluvial_plot_first_state <- 
  create_sankey_first_state(trajectories_limit = trajectories_limit_counts)

alluvial_plot_first_state

# # Save the plot
# ggsave(
#   plot = alluvial_plot_first_state,
#   filename = paste0(wd, "Outputs/", "overall_switches_alluvial_first_state", 
#                     today(), ".png"),
#   dpi = 300, width = 13, height = 8, bg = "white"
# )


#===================== 60-day interval ========================================= 

prescriptions_combs_final_60 <- get_treatment_intervals(
  prescriptions_prep = prescriptions_prep, overlap_days = 60, 
  interrupt_cutoff_intervals = interrupt_cutoff_intervals)

trajectories_limit_60 <- derive_trajectories(
  prescriptions_combs_final = prescriptions_combs_final_60, 
  max_switches = max_switches)

# n = 36705 after suppressing small counts
trajectories_limit_counts_60 <- make_safe_counts(trajectories_limit_60)

alluvial_plot_60 <- create_sankey(trajectories_limit = trajectories_limit_counts_60)

alluvial_plot_60

# # Save the plot
# ggsave(
#   plot = alluvial_plot_60,
#   filename = paste0(wd, "Outputs/", "overall_switches_alluvial_60day", today(), ".png"),
#   dpi = 300, width = 13, height = 9, bg = "white"
# )

# # Save as PDF
# ggsave(
#   plot = alluvial_plot_60,
#   filename = paste0(wd, "Outputs/", "overall_switches_alluvial_60day", today(), ".pdf"),
#   dpi = 300, width = 13, height = 9, bg = "white"
# )

#===================== 1-day interval w/ 90-days as interruption ===============

prescriptions_combs_final_1 <- get_treatment_intervals(
  prescriptions_prep = prescriptions_prep, overlap_days = 1, 
  interrupt_cutoff_intervals = 90)

trajectories_limit_1 <- derive_trajectories(
  prescriptions_combs_final = prescriptions_combs_final_1, 
  max_switches = max_switches)

# n = 37085 after suppressing small counts
trajectories_limit_counts_1 <- make_safe_counts(trajectories_limit_1)

alluvial_plot_1 <- create_sankey(trajectories_limit = trajectories_limit_counts_1)

alluvial_plot_1

# # Save the plot
# ggsave(
#   plot = alluvial_plot_1,
#   filename = paste0(wd, "Outputs/", "overall_switches_alluvial_1day", today(), ".png"),
#   dpi = 300, width = 13, height = 9, bg = "white"
# )

# # Save as PDF
# ggsave(
#   plot = alluvial_plot_1,
#   filename = paste0(wd, "Outputs/", "overall_switches_alluvial_1day", today(), ".pdf"),
#   dpi = 300, width = 13, height = 9, bg = "white"
# )

#===================== Lithium focus ===========================================

#Lithium category

prescriptions_prep_li <- prescriptions_prep  %>%
  mutate(drug_type = if_else(grepl("Lithium", ingredient, ignore.case = TRUE),
                             "Lithium",
                             drug_type))

prescriptions_combs_final_li <- get_treatment_intervals(
  prescriptions_prep = prescriptions_prep_li, overlap_days = 30, 
  interrupt_cutoff_intervals = interrupt_cutoff_intervals)

# Collapse lithium + 2 or more classes into one category
prescriptions_combs_final_li2 <- prescriptions_combs_final_li %>%
  mutate(
    # Preserve the raw combo only for non-lithium rows
    nonli_combo_raw = if_else(!grepl("Lithium", drug_combination, ignore.case = TRUE),
                              drug_combination, NA_character_),
    drug_combination = case_when(
      drug_combination %in% c("AD; AP; Lithium", "AD; AP; Lithium; MS",
                              "AD; Lithium; MS", "AP; Lithium; MS") ~ "Lithium; 2+",
      drug_combination %in% c("AD", "AD; AP", "AD; AP; MS", "AD; MS", "AP", 
                              "AP; MS", "MS") ~ "Non-Lithium",
      # Reorder some names so lithium comes first
      drug_combination == "AD; Lithium" ~ "Lithium; AD",
      drug_combination == "AP; Lithium" ~ "Lithium; AP",
      .default = drug_combination))

table(prescriptions_combs_final_li2$drug_combination)

trajectories_limit_li <- derive_trajectories(
  prescriptions_combs_final = prescriptions_combs_final_li2, 
  max_switches = max_switches, 
  lithium = TRUE)

# # Save for Shiny app
# n = 36955 after suppressing small counts
trajectories_limit_li_counts <- make_safe_counts(trajectories_limit_li)
# saveRDS(trajectories_limit_li_counts, 
#         paste0("bpd-prescribing-uk-shiny/data/",
#                "trajectories_limit_li_counts.rds"))


## Define factor order and colors
# Order of states in each step
state_levels_li <- c(
  "Lithium", "Lithium; AD", "Lithium; AP", "Lithium; MS", "Lithium; 2+", 
  "Non-Lithium", "Switch Non-Lithium", "Int/Disc", 
  "No Switch"
)

# Define colors
state_colors_li <- c(
  "Lithium" = "#FFC33B",          # bright red
  "Lithium; AD" = "#E20134",      # purple
  "Lithium; AP" = "#009f81",      # orange
  "Lithium; MS" = "#FF6E3A",      # green
  "Lithium; 2+" = "#994F00",  # brown
  "Non-Lithium" = "#008DF9",          # yellow
  "Switch Non-Lithium" = "#00C2F9",
  "Int/Disc" = "#999999",    # gray
  "No Switch" = "#F5F5F5"
)

# Map abbreviation states
state_full_labels_li <- c(
  "Lithium" = "Lithium", 
  "Lithium; AD" = "Lithium & Antidepressants", 
  "Lithium; AP" = "Lithium & Antipsychotics", 
  "Lithium; MS" = "Lithium & Mood Stabilisers",
  "Lithium; 2+" = "Lithium & 2 or More Classes",
  "Non-Lithium" = "Non-Lithium",
  "Switch Non-Lithium" = "Switch Within Non-Lithium",
  "Int/Disc" = "Interruption/Discontinuation",    
  "No Switch" = "No Switch")

alluvial_plot_li <- create_sankey(trajectories_limit = trajectories_limit_li_counts,
                                  state_levels = state_levels_li, 
                                  state_colors = state_colors_li,
                                  state_full_labels = state_full_labels_li, 
                                  text_size = 2.5)

alluvial_plot_li

# # Save the plot
# ggsave(
#   plot = alluvial_plot_li,
#   filename = paste0(wd, "Outputs/", "overall_switches_alluvial_lithium", 
#                     today(), ".png"),
#   dpi = 300, width = 13, height = 9, bg = "white"
# )

# # Save as PDF
# ggsave(
#   plot = alluvial_plot_li,
#   filename = paste0(wd, "Outputs/", "overall_switches_alluvial_lithium", 
#                     today(), ".pdf"),
#   dpi = 300, width = 13, height = 9, bg = "white"
# )

# Lithium colored by first state

alluvial_plot_li_first_state <- create_sankey_first_state(
  trajectories_limit = trajectories_limit_li_counts,
  state_levels = state_levels_li, 
  state_colors = state_colors_li,
  state_full_labels = state_full_labels_li, alpha_input = 0.8)

alluvial_plot_li_first_state

# # Save the plot
# ggsave(
#   plot = alluvial_plot_li_first_state,
#   filename = paste0(wd, "Outputs/", "overall_switches_alluvial_lithium_first_state",
#                     today(), ".png"),
#   dpi = 300, width = 13, height = 9, bg = "white"
# )


# Lithium with all categories

#Lithium category

prescriptions_prep_li_all <- prescriptions_prep  %>%
  mutate(drug_type = if_else(grepl("Lithium", ingredient, ignore.case = TRUE),
                             "Lithium",
                             drug_type))

prescriptions_combs_final_li_all <- get_treatment_intervals(
  prescriptions_prep = prescriptions_prep_li_all, overlap_days = 30, 
  interrupt_cutoff_intervals = interrupt_cutoff_intervals)

# Collapse no lithium into one category
prescriptions_combs_final_li_all2 <- prescriptions_combs_final_li_all %>%
  mutate(drug_combination = case_when(
    drug_combination %in% c("AD", "AD; AP", "AD; AP; MS", "AD; MS", "AP", 
                            "AP; MS", "MS") ~ "Non-Lithium",
    # Reorder some names so lithium comes first
    drug_combination == "AD; Lithium" ~ "Lithium; AD",
    drug_combination == "AP; Lithium" ~ "Lithium; AP",
    drug_combination == "AD; AP; Lithium" ~ "Lithium; AD; AP",
    drug_combination == "AD; Lithium; MS" ~ "Lithium; AD; MS",
    drug_combination == "AP; Lithium; MS" ~ "Lithium; AP; MS",
    drug_combination == "AD; AP; Lithium; MS" ~ "Lithium; AD; AP; MS",
    .default = drug_combination
  ))

table(prescriptions_combs_final_li_all2$drug_combination)

trajectories_limit_li_all <- derive_trajectories(
  prescriptions_combs_final = prescriptions_combs_final_li_all2, 
  max_switches = max_switches)

## Define factor order and colors
# Order of states in each step
state_levels_li_all <- c(
  "Lithium", "Lithium; AD", "Lithium; AP", "Lithium; MS", 
  "Lithium; AD; AP", "Lithium; AD; MS", "Lithium; AP; MS", 
  "Lithium; AD; AP; MS",
  "Non-Lithium", "Int/Disc", 
  "No Switch"
)

# Define colors
state_colors_li_all <- c(
  "Lithium" = "#D35FB7",          # bright red
  "Lithium; AD" = "#E41A1C",      # purple
  "Lithium; AP" = "#377EB8",      # orange
  "Lithium; MS" = "#FFD92F",      # green
  "Lithium; AD; AP" = "#984EA3", 
  "Lithium; AD; MS" = "#FF7F00", 
  "Lithium; AP; MS" = "#4DAF4A",
  "Lithium; AD; AP; MS" = "#994F00",
  "Non-Lithium" = "#40B0A6",         
  "Int/Disc" = "#999999",    # gray
  "No Switch" = "#F5F5F5"
)

# Map abbreviation states
state_full_labels_li_all <- c(
  "Lithium" = "Lithium", 
  "Lithium; AD" = "Lithium & Antidepressants", 
  "Lithium; AP" = "Lithium & Antipsychotics", 
  "Lithium; MS" = "Lithium & Mood Stabilisers",
  "Lithium; AD; AP" = "Lithium & Antidepressants & Antipsychotics",
  "Lithium; AD; MS" = "Lithium & Antidepressants & Mood Stabilisers",
  "Lithium; AP; MS" = "Lithium & Antipsychotics & Mood Stabilisers",
  "Lithium; AD; AP; MS" = "Lithium & Antidepressants & \nAntipsychotics & Mood Stabilisers",
  "Non-Lithium" = "Non-Lithium",
  "Int/Disc" = "Interruption/Discontinuation",    
  "No Switch" = "No Switch")

alluvial_plot_li_all <- create_sankey(trajectories_limit = trajectories_limit_li_all,
                                      state_levels = state_levels_li_all, 
                                      state_colors = state_colors_li_all,
                                      state_full_labels = state_full_labels_li_all)

alluvial_plot_li_all


#===================== Antidepressant first stage only =========================

trajectories_limit_AD <- trajectories_limit_counts %>%
  filter(state_1 == "AD")
# filter(grepl("AD", state_1))

alluvial_plot_AD <- create_sankey(trajectories_limit = trajectories_limit_AD)

alluvial_plot_AD

# # Save the plot
# ggsave(
#   plot = alluvial_plot_AD,
#   filename = paste0(wd, "Outputs/", "overall_switches_alluvial_AD", today(), ".png"),
#   dpi = 300, width = 13, height = 9, bg = "white"
# )


#================== Heatmap initial and second line ============================

# Calculate transition percentages
# n = 36910
trajectories_limit_new <- trajectories_limit_counts %>%
  mutate(state_2 = case_when(is.na(state_2) ~ "No Switch", TRUE ~ state_2))

# Get counts, summing over state_3
if (!"n" %in% names(trajectories_limit_new)) {
  trajectories_limit_new <- count(trajectories_limit_new, state_1, state_2, 
                                  name = "n")
} else {
  trajectories_limit_new <- trajectories_limit_new %>% 
    group_by(., state_1, state_2) %>%
    summarise(n = sum(n, na.rm = TRUE), .groups = "drop")
}
transition_matrix <- trajectories_limit_new %>%
  group_by(state_1) %>%
  mutate(percentage = (n / sum(n)) * 100, total_n = sum(n)) %>%
  ungroup()

# Specify exact order for both axes
axis_order <- c("AD", "AP", "MS", "AD; AP", "AD; MS", "AP; MS", 
                "AD; AP; MS", "No Switch", "Int/Disc")

# Calculate grand total and create labels with percentages
grand_total <- sum(unique(transition_matrix$total_n))
sample_sizes <- transition_matrix %>%
  distinct(state_1, total_n) %>%
  mutate(
    state_1 = factor(state_1, levels = axis_order),
    pct_of_total = (total_n / grand_total) * 100,
    label_simple = paste0(state_1, "\nn=", format(total_n, big.mark = ","), 
                          " (", sprintf("%.1f%%", pct_of_total), ")")) %>%
  arrange(state_1)

transition_matrix_final <- transition_matrix %>%
  left_join(sample_sizes %>% select(state_1, label_simple), by = "state_1") %>%
  mutate(
    state_2 = factor(state_2, levels = axis_order),
    label_simple = factor(label_simple, levels = sample_sizes$label_simple))

# Calculate dynamic max for color scale
max_percentage <- max(transition_matrix_final$percentage, na.rm = TRUE)
color_max <- ceiling(max_percentage / 10) * 10 + 10

# Base theme for both plots
base_theme <- theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 8, lineheight = 0.9),
    axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 10)),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 10)))

# Plot 1: Fixed row heights
initial_transition <- ggplot(transition_matrix_final, aes(x = state_2, y = label_simple, fill = percentage)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            color = "black", size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#f7fbff", mid = "#4292c6", high = "#08519c",
                       midpoint = color_max / 2, limits = c(0, color_max),
                       name = "Percentage (%)") +
  labs(x = "Second-Line Treatment", y = "First-Line Treatment"
       #   title = "Treatment Transitions (Fixed Row Heights)"
  ) +
  base_theme +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = 15, barheight = 0.8))

# Plot 2: Proportional row heights
sample_sizes_positions <- sample_sizes %>%
  mutate(
    plot_height = pmax(pct_of_total, 3),
    y_max = cumsum(plot_height),
    y_min = lag(y_max, default = 0),
    y_center = (y_max + y_min) / 2)

transition_matrix_proportional <- transition_matrix_final %>%
  left_join(sample_sizes_positions %>% select(state_1, y_min, y_max, y_center), by = "state_1") %>%
  mutate(
    x_num = as.numeric(state_2),
    x_min = x_num - 0.45,
    x_max = x_num + 0.45)

initial_transition_proportional <- ggplot(transition_matrix_proportional) +
  geom_rect(aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max, fill = percentage),
            color = "white", linewidth = 0.8) +
  geom_text(aes(x = x_num, y = y_center, label = sprintf("%.1f%%", percentage)), 
            color = "black", size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#f7fbff", mid = "#4292c6", high = "#08519c",
                       midpoint = color_max / 2, limits = c(0, color_max),
                       name = "Percentage (%)") +
  scale_x_continuous(breaks = 1:length(axis_order), labels = axis_order, expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = sample_sizes_positions$y_center, 
                     labels = sample_sizes_positions$label_simple, expand = c(0, 0)) +
  labs(x = "Second-Line Treatment", y = "First-Line Treatment"
       # title = "Treatment Transitions (Proportional Row Heights)"
  ) +
  base_theme +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = 15, barheight = 0.8))

# Display both plots
print(initial_transition)
print(initial_transition_proportional)


# ================== Heatmap for Lithium =======================================

# Calculate transition percentages
trajectories_limit_new_li <- trajectories_limit_li_counts %>%
  mutate(state_2 = case_when(is.na(state_2) ~ "No Switch", TRUE ~ state_2))

# Get counts, summing over state_3
if (!"n" %in% names(trajectories_limit_new_li)) {
  trajectories_limit_new_li <- count(trajectories_limit_new_li, state_1, state_2, 
                                  name = "n")
} else {
  trajectories_limit_new_li <- trajectories_limit_new_li %>% 
    group_by(., state_1, state_2) %>%
    summarise(n = sum(n, na.rm = TRUE), .groups = "drop")
}
transition_matrix_li <- trajectories_limit_new_li %>%
  group_by(state_1) %>%
  mutate(percentage = (n / sum(n)) * 100, total_n = sum(n)) %>%
  ungroup()

# Specify exact order for both axes
axis_order_li <- c(
  "Lithium", "Lithium; AD", "Lithium; AP", "Lithium; MS", "Lithium; 2+", 
  "Non-Lithium", "Switch Non-Lithium", "No Switch", "Int/Disc")

# Calculate grand total and create labels with percentages
grand_total_li <- sum(unique(transition_matrix_li$total_n))
sample_sizes_li <- transition_matrix_li %>%
  distinct(state_1, total_n) %>%
  mutate(
    state_1 = factor(state_1, levels = axis_order_li),
    pct_of_total = (total_n / grand_total_li) * 100,
    label_simple = paste0(state_1, "\nn=", format(total_n, big.mark = ","), 
                          " (", sprintf("%.1f%%", pct_of_total), ")")) %>%
  arrange(state_1)

transition_matrix_final_li <- transition_matrix_li %>%
  left_join(sample_sizes_li %>% select(state_1, label_simple), by = "state_1") %>%
  mutate(
    state_2 = factor(state_2, levels = axis_order_li),
    label_simple = factor(label_simple, levels = sample_sizes_li$label_simple))

# Calculate dynamic max for color scale
max_percentage_li <- max(transition_matrix_final_li$percentage, na.rm = TRUE)
color_max_li <- ceiling(max_percentage_li / 10) * 10 + 10

# Plot 1: Fixed row heights (Lithium)
initial_transition_li <- ggplot(transition_matrix_final_li, aes(x = state_2, y = label_simple, fill = percentage)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            color = "black", size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#f7fbff", mid = "#4292c6", high = "#08519c",
                       midpoint = color_max_li / 2, limits = c(0, color_max_li),
                       name = "Percentage (%)") +
  labs(x = "Second-Line Treatment", y = "First-Line Treatment"
       #   title = "Lithium Treatment Transitions (Fixed Row Heights)"
  ) +
  base_theme +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = 15, barheight = 0.8))

# Plot 2: Proportional row heights (Lithium)
sample_sizes_positions_li <- sample_sizes_li %>%
  mutate(
    plot_height = pmax(pct_of_total, 3),
    y_max = cumsum(plot_height),
    y_min = lag(y_max, default = 0),
    y_center = (y_max + y_min) / 2)

transition_matrix_proportional_li <- transition_matrix_final_li %>%
  left_join(sample_sizes_positions_li %>% select(state_1, y_min, y_max, y_center), by = "state_1") %>%
  mutate(
    x_num = as.numeric(state_2),
    x_min = x_num - 0.45,
    x_max = x_num + 0.45)

initial_transition_proportional_li <- ggplot(transition_matrix_proportional_li) +
  geom_rect(aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max, fill = percentage),
            color = "white", linewidth = 0.8) +
  geom_text(aes(x = x_num, y = y_center, label = sprintf("%.1f%%", percentage)), 
            color = "black", size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#f7fbff", mid = "#4292c6", high = "#08519c",
                       midpoint = color_max_li / 2, limits = c(0, color_max_li),
                       name = "Percentage (%)") +
  scale_x_continuous(breaks = 1:length(axis_order_li), labels = axis_order_li, expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = sample_sizes_positions_li$y_center, 
                     labels = sample_sizes_positions_li$label_simple, expand = c(0, 0)) +
  labs(x = "Second-Line Treatment", y = "First-Line Treatment"
       # title = "Lithium Treatment Transitions (Proportional Row Heights)"
  ) +
  base_theme +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = 15, barheight = 0.8))

# Display lithium plots
print(initial_transition_li)
print(initial_transition_proportional_li)

# Save heatmaps
plots <- list(
  initial_transition = initial_transition,
  initial_transition_proportional = initial_transition_proportional,
  initial_transition_li = initial_transition_li,
  initial_transition_proportional_li = initial_transition_proportional_li)

for (nm in names(plots)) {
  ggsave(
    filename = paste0(nm, "_heatmap_", today(), ".pdf"),
    path = file.path(wd, "Outputs"),
    plot = plots[[nm]],
    width = 8, height = 5, units = "in")}

# # Save some as png
# ggsave(
#   plot = initial_transition,
#   filename = paste0(wd, "Outputs/", "initial_transition_heatmap_", today(), ".png"),
#   dpi = 300, width = 8, height = 5, units = "in", bg = "white"
# )
# ggsave(
#   plot = initial_transition_li,
#   filename = paste0(wd, "Outputs/", "initial_transition_li_heatmap_", today(), ".png"),
#   dpi = 300, width = 8, height = 5, units = "in", bg = "white"
# )

# # Same some as PDF
# ggsave(
#   plot = initial_transition,
#   filename = paste0(wd, "Outputs/", "initial_transition_heatmap_", today(), ".pdf"),
#   dpi = 300, width = 8, height = 5, units = "in", bg = "white"
# )
# ggsave(
#   plot = initial_transition_li,
#   filename = paste0(wd, "Outputs/", "initial_transition_li_heatmap_", today(), ".pdf"),
#   dpi = 300, width = 8, height = 5, units = "in", bg = "white"
# )
rm(plots)



#================== Heatmap second and third line ==============================

# -------- Second → Third transitions (overall) --------

# Ensure explicit labels for missing states
trajectories_limit_23 <- trajectories_limit_counts %>%
  mutate(
    state_2 = ifelse(is.na(state_2), "No Switch", state_2),
    state_3 = ifelse(is.na(state_3), "No Switch", state_3)
  )

# Aggregate second → third, then compute within–state_2 percentages
if (!"n" %in% names(trajectories_limit_23)) {
  trajectories_limit_23 <- trajectories_limit_23 %>%
    group_by(state_2, state_3) %>%
    summarise(n = n(), .groups = "drop")
} else {
  trajectories_limit_23 <- trajectories_limit_23 %>% 
    group_by(state_2, state_3) %>%
    summarise(n = sum(n, na.rm = TRUE), .groups = "drop")
}
transition_matrix_23 <- trajectories_limit_23 %>%
  group_by(state_2) %>%
  mutate(percentage = (n / sum(n)) * 100, total_n = sum(n)) %>%
  ungroup()

# Order for both axes
axis_order <- c("AD", "AP", "MS", "AD; AP", "AD; MS", "AP; MS",
                "AD; AP; MS", "No Switch", "Int/Disc")

# Row labels: include n and % of grand total (based on second-line totals)
grand_total_23 <- sum(unique(transition_matrix_23$total_n), na.rm = TRUE)

sample_sizes_23 <- transition_matrix_23 %>%
  distinct(state_2, total_n) %>%
  mutate(
    state_2      = factor(state_2, levels = axis_order),
    pct_of_total = if (grand_total_23 > 0) (total_n / grand_total_23) * 100 else NA_real_,
    label_simple = paste0(state_2, "\nn=", format(total_n, big.mark = ","), 
                          " (", sprintf("%.1f%%", pct_of_total), ")")
  ) %>%
  arrange(state_2)

transition_matrix_final_23 <- transition_matrix_23 %>%
  left_join(sample_sizes_23 %>% select(state_2, label_simple), by = "state_2") %>%
  mutate(
    state_3     = factor(state_3, levels = axis_order),
    label_simple = factor(label_simple, levels = sample_sizes_23$label_simple)
  )

# Hide "No Switch" row only at plot time (keep percentages unchanged) 
transition_matrix_plot <- transition_matrix_final_23 %>%
  dplyr::filter(state_2 != "No Switch") %>%
  droplevels()  # drop the now-empty y-factor level

# Dynamic max for color scale
max_percentage_23 <- max(transition_matrix_plot$percentage, na.rm = TRUE)
color_max_23 <- ceiling(max_percentage_23 / 5) * 5 + 5

# Plot: Fixed row heights (Second → Third)
second_to_third <- ggplot(transition_matrix_plot,
                          aes(x = state_3, y = label_simple, fill = percentage)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)),
            color = "black", size = 3.5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#f7fbff", mid = "#4292c6", high = "#08519c",
    midpoint = color_max_23 / 2, limits = c(0, color_max_23),
    name = "Percentage (%)"
  ) +
  labs(
    x = "Third-Line Treatment",
    y = "Second-Line Treatment"
    # title = "Second → Third (Fixed Row Heights)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 10),
    axis.text.y  = element_text(face = "bold", size = 8, lineheight = 0.9),
    axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 10)),
    panel.grid   = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8),
    plot.title   = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 10))
  ) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5,
                               barwidth = 15, barheight = 0.8))

second_to_third

# Save output
# ggsave(
#   plot = second_to_third,
#   filename = paste0(wd, "Outputs/", "second_transition_heatmap_", today(), ".png"),
#   dpi = 300, width = 8, height = 5, units = "in", bg = "white"
# )

# Save as PDF
# ggsave(
#   plot = second_to_third,
#   filename = paste0(wd, "Outputs/", "second_transition_heatmap_", today(), ".pdf"),
#   dpi = 300, width = 8, height = 5, units = "in", bg = "white"
# )

# -------- Second → Third transitions (Lithium-focused) --------

trajectories_limit_li_23 <- trajectories_limit_li_counts %>%
  mutate(
    state_2 = ifelse(is.na(state_2), "No Switch", state_2),
    state_3 = ifelse(is.na(state_3), "No Switch", state_3)
  )
# Aggregate second → third, then compute within–state_2 percentages
if (!"n" %in% names(trajectories_limit_li_23)) {
  trajectories_limit_li_23 <- trajectories_limit_li_23 %>%
    group_by(state_2, state_3) %>%
    summarise(n = n(), .groups = "drop")
} else {
  trajectories_limit_li_23 <- trajectories_limit_li_23 %>% 
    group_by(state_2, state_3) %>%
    summarise(n = sum(n, na.rm = TRUE), .groups = "drop")
}
transition_matrix_23_li <- trajectories_limit_li_23 %>%
  group_by(state_2) %>%
  mutate(percentage = (n / sum(n)) * 100, total_n = sum(n)) %>%
  ungroup()

axis_order_li <- c(
  "Lithium", "Lithium; AD", "Lithium; AP", "Lithium; MS", "Lithium; 2+",
  "Non-Lithium", "Switch Non-Lithium", "No Switch", "Int/Disc"
)

grand_total_23_li <- sum(unique(transition_matrix_23_li$total_n), na.rm = TRUE)

sample_sizes_23_li <- transition_matrix_23_li %>%
  distinct(state_2, total_n) %>%
  mutate(
    state_2      = factor(state_2, levels = axis_order_li),
    pct_of_total = if (grand_total_23_li > 0) (total_n / grand_total_23_li) * 100 else NA_real_,
    label_simple = paste0(state_2, "\nn=", format(total_n, big.mark = ","), 
                          " (", sprintf("%.1f%%", pct_of_total), ")")
  ) %>%
  arrange(state_2)

transition_matrix_final_23_li <- transition_matrix_23_li %>%
  left_join(sample_sizes_23_li %>% select(state_2, label_simple), by = "state_2") %>%
  mutate(
    state_3      = factor(state_3, levels = axis_order_li),
    label_simple = factor(label_simple, levels = sample_sizes_23_li$label_simple)
  )

# Hide "No Switch" row only at plot time (keep percentages unchanged) 
transition_matrix_plot_li <- transition_matrix_final_23_li %>%
  dplyr::filter(state_2 != "No Switch") %>%
  droplevels()  # drop the now-empty y-factor level

max_percentage_23_li <- max(transition_matrix_plot_li$percentage, na.rm = TRUE)
color_max_23_li <- ceiling(max_percentage_23_li / 5) * 5 + 5

second_to_third_li <- ggplot(transition_matrix_plot_li,
                             aes(x = state_3, y = label_simple, fill = percentage)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)),
            color = "black", size = 3.5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#f7fbff", mid = "#4292c6", high = "#08519c",
    midpoint = color_max_23_li / 2, limits = c(0, color_max_23_li),
    name = "Percentage (%)"
  ) +
  labs(
    x = "Third-Line Treatment",
    y = "Second-Line Treatment"
    # title = "Second → Third (Lithium, Fixed Row Heights)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 10),
    axis.text.y  = element_text(face = "bold", size = 8, lineheight = 0.9),
    axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 10)),
    panel.grid   = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8),
    plot.title   = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 10))
  ) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5,
                               barwidth = 15, barheight = 0.8))

second_to_third_li

# Save output
# ggsave(
#   plot = second_to_third_li,
#   filename = paste0(wd, "Outputs/", "second_transition_li_heatmap_", today(), ".png"),
#   dpi = 300, width = 8, height = 5, units = "in", bg = "white"
# )

# Save as PDF
# ggsave(
#   plot = second_to_third_li,
#   filename = paste0(wd, "Outputs/", "second_transition_li_heatmap_", today(), ".pdf"),
#   dpi = 300, width = 8, height = 5, units = "in", bg = "white"
# )
