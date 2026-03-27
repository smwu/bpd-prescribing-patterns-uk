#===============================================================================
# Create line plots of the proportion of patients prescribed each medication 
# class in the 12 months before and after first-recorded BPD diagnosis
# 
# Data:
#   - Data/bipolar_cohort_ads.rdata: Antidepressant prescriptions
#   - Data/bipolar_cohort_moodstabs.rdata: Mood stabiliser prescriptions
#   - Data/bipolar_cohort_aps.rdata: Antipsychotic prescriptions
#   - Data/bipolar_cohort.rdata: Patient information 
# 
# Outputs:
#   - Outputs/Percentage_lineplot: Line plot of monthly medication prescription percentages  
#   - Outputs/Percentage_lineplot_pre_covid: Prescription percentages restricted to pre-Covid
#   - Outputs/Percentage_lineplot_from_covid: Prescription percentages restricted to post-Covid
#   - Outputs/Percentage_lineplot_pre_from_covid_combined: Combined pre-Covid and post-Covid plots
# ==============================================================================

# Clear memory
rm(list = ls())

library(lubridate)
library(gtsummary)
library(tidyverse)
library(ggforce)
library(gt)
library(ggrepel)
library(ggpubr)
library(rlang)
library(flextable)
library(officer)

load("Data/bipolar_cohort_ads.rdata")
load("Data/bipolar_cohort_moodstabs.rdata")
load("Data/bipolar_cohort_aps.rdata")
load("Data/bipolar_cohort.rdata")

wd <- getwd()

#======= Prescription data

# Combine valproate variants into "Valproate"
bipolar_cohort_moodstabs <- bipolar_cohort_moodstabs %>%
  mutate(ingredient = case_when(
    ingredient %in% c("Valproate", "Valproate/Valproic", "Valproic") ~ "Valproate",
    TRUE ~ ingredient
  ))

# Get first bipolar diagnosis date
bipolar_diagnosis_date <- bipolar_cohort %>%
  mutate(bipolar_diag_date = as.Date(first_bipolar_date)) %>%
  select(patid, bipolar_diag_date)

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

# Join with diagnosis data: 40,033 patients appear in both diagnosis and 
# prescription dataset
prescriptions_with_diag <- all_prescriptions %>%
  left_join(bipolar_diagnosis_date, by = "patid")

# Filter by time window (within 1 year before or after diagnosis)
prescriptions_filtered_oneyr_beforeafter <- prescriptions_with_diag %>%
  mutate(diff_days = as.numeric(issuedate - bipolar_diag_date)) %>%
  filter(abs(diff_days) <= 365)

# Number of patients with prescription data within one year of diagnosis: 38,446
# So 2519 with BP do not have any prescription data within one year
length(unique(prescriptions_filtered_oneyr_beforeafter$patid))

# Number of patients with prescription data within one year *before* diagnosis: 32,460
prescriptions_filtered_oneyr_before <- prescriptions_filtered_oneyr_beforeafter %>% 
  filter(diff_days < 0)
length(unique(prescriptions_filtered_oneyr_before$patid))

# Number of patients with prescription data within one year *after* diagnosis: 37,208
# So 3757 with BP do not have any prescription data within one year after diagnosis
prescriptions_filtered_oneyr_after <- prescriptions_filtered_oneyr_beforeafter %>% 
  filter(diff_days >= 0)
length(unique(prescriptions_filtered_oneyr_after$patid))


## Calculate months since diagnosis (-12:-1, 1:12)
# Date of diagnosis is the baseline for the calculation
# Month intervals are by 30.42
# After diagnosis: 0 (inclusive) to 30.42 days after is considered month 1, 
# then increases by 30.42 days up to month 12 (334.63 to 365.04) after
# Before diagnosis: -30.42 to 0 (exclusive) is considered month -1, then 
# up to month -12 (-334.63 to -365.04)
prescriptions_with_month <- prescriptions_filtered_oneyr_beforeafter %>%
  mutate(
    issuedate = as.Date(issuedate),
    bipolar_diag_date = as.Date(bipolar_diag_date),
    month_since_diag = case_when(
      diff_days == 0 ~ 1,  # On date of diagnosis => consider as month 1
      diff_days > 0 ~ ceiling(diff_days / 30.42), # After diagnosis (month 1, 2, 3,...)
      diff_days < 0 ~ floor(diff_days / 30.42), # Before diagnosis (month -1, -2, -3,...)
      TRUE ~ NA_real_
    )) %>%
  filter(month_since_diag >= -12 & month_since_diag <= 12)

# Test month calculation
ceiling(c(1, 29, 30, 31, 364, 365, 334.63, 365.04) / 30.42)
floor(-c(1, 29, 30, 31, 364, 365, 334.63, 365.04) / 30.42)

# Convert month_since_drug to factor
prescriptions_with_month <- prescriptions_with_month %>% 
  mutate(month_since_diag = factor(
    month_since_diag,
    levels = c(-12:-1, 1:12),            
    ordered = TRUE
  )) 


## Drug Type ##
# For each patient and each month since diagnosis, combine all unique drug 
# classes prescribed into a single semicolon-separated list
monthly_drugs <- prescriptions_with_month %>%
  group_by(patid, month_since_diag) %>%
  summarise(
    drug_combined = paste(unique(drug_type), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(patid, month_since_diag) 

# Pivot to wide format — one row per patid
wide_monthly_drugs <- monthly_drugs %>%
  pivot_wider(
    names_from = month_since_diag,
    values_from = drug_combined,
    names_prefix = "month_"
  )

# Add missing month columns (month_-12 to month_12) and replace NAs with "None"
all_months <- paste0("month_", c(-12:-1, 1:12))
wide_monthly_drugs <- wide_monthly_drugs %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(all_of(setdiff(all_months, names(.))), ~NA_character_)) %>%
  mutate(across(all_of(all_months), ~replace_na(.x, "None")))


#========= Incorporate prescription data into full cohort

# Ensure every unique patid is included (even if no prescriptions at all)
unique_patids <- bipolar_cohort %>%
  select(patid) %>%
  distinct()

final_df_monthly_drugs <- unique_patids %>%
  left_join(wide_monthly_drugs, by = "patid") %>%
  mutate(across(starts_with("month_"), ~replace_na(.x, "None")))

# Reorder columns
final_df_monthly_drugs <- final_df_monthly_drugs %>%
  select(patid, all_of(all_months))


#=========== Bipolar cohort data censoring

# Censor at earliest of death, lcd, or regendate
censor_dates <- bipolar_cohort %>%
  select(patid, deathdate, regenddate, lcd) %>%
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
  select(patid, earliest_exit, diff_since_diag, month_since_diag, within_12_months)

# Incorporate censoring information into final dataset
# All months after censoring are set to NA
# Note: "None" means no prescription recorded, whereas NA means censored
final_df <- final_df_monthly_drugs %>%
  left_join(
    censor_dates %>% select(patid, month_since_diag, within_12_months),
    by = "patid"
  ) %>%
  mutate(across(all_of(paste0("month_", 1:12)), ~as.character(.x))) %>%
  rowwise() %>%
  mutate(across(
    paste0("month_", 1:12),
    ~if (!is.na(month_since_diag) && within_12_months &&
         as.integer(str_extract(cur_column(), "\\d+")) >= month_since_diag) {
      NA_character_
    } else {
      .x
    }
  )) %>%
  ungroup() %>%
  select(-month_since_diag, -within_12_months)

# Check final data
table(final_df$`month_-1`, useNA = "always") # No NAs before diagnosis
table(final_df$`month_1`, useNA = "always")


# Summarize how many meeting the criteria for censoring
censor_summary <- bipolar_cohort %>%
  select(patid, deathdate, regenddate, lcd) %>%
  left_join(bipolar_diagnosis_date, by = "patid") %>%
  mutate(
    deathdate = as.Date(deathdate),
    regenddate = as.Date(regenddate),
    lcd = as.Date(lcd),
    bipolar_diag_date = as.Date(bipolar_diag_date),
    death_within_12mo = !is.na(deathdate) & deathdate >= bipolar_diag_date & 
      deathdate <= bipolar_diag_date + 365,
    regend_within_12mo = !is.na(regenddate) & regenddate >= bipolar_diag_date & 
      regenddate < bipolar_diag_date + 365,
    lcd_within_12mo = !is.na(lcd) & lcd >= bipolar_diag_date & 
      lcd < bipolar_diag_date + 365,
    any_within_12mo = death_within_12mo | regend_within_12mo | lcd_within_12mo
  )

censor_counts <- censor_summary %>%
  summarise(
    death_within_12mo_n = sum(death_within_12mo, na.rm = TRUE),
    regend_within_12mo_n = sum(regend_within_12mo, na.rm = TRUE),
    lcd_within_12mo_n = sum(lcd_within_12mo, na.rm = TRUE),
    total_patients_within_12mo = sum(any_within_12mo, na.rm = TRUE)
  )

print(censor_counts)

# Check if there were any censored due to lcd or regend rather than death
temp <- censor_summary %>% 
  filter(any_within_12mo == TRUE) %>% 
  filter(!(death_within_12mo == TRUE))


# ========== Overall % of patients per drug combo pre- and post-diagnosis ======

# Long format with month number and period (pre vs post)
long_year <- final_df %>%
  pivot_longer(
    cols = starts_with("month_"),
    names_to = "month",
    values_to = "drug_combination"
  ) %>%
  mutate(
    month_num = as.integer(str_extract(month, "-?\\d+")),
    period = case_when(
      month_num < 0 ~ "Year pre-diagnosis",
      month_num > 0 ~ "Year post-diagnosis",
      TRUE ~ NA_character_   # should drop the 0 month if it ever existed
    )
  ) %>%
  # Drop censored months (NA) and months outside the ±12 window (already handled)
  filter(!is.na(period), !is.na(drug_combination)) %>%
  # Recode "None" to match your plotting labels
  mutate(
    drug_combination = factor(
      drug_combination,
      levels = c("None", "AD", "AP", "MS",
                 "AD; AP", "AD; MS", "AP; MS", "AD; AP; MS"),
      labels = c("No AD, AP, or MS", "AD", "AP", "MS",
                 "AD; AP", "AD; MS", "AP; MS", "AD; AP; MS")
    )
  )

# Denominator: number of patients with any observed time in each period
denom_period <- long_year %>%
  distinct(patid, period) %>%
  count(period, name = "n_patients_period")

# For each period, count patients who ever had each drug combination at least once
overall_yearly_pct <- long_year %>%
  distinct(patid, period, drug_combination) %>%   # "ever on" that combo in that year
  count(period, drug_combination, name = "n_patients") %>%
  left_join(denom_period, by = "period") %>%
  mutate(
    percent_patients = 100 * n_patients / n_patients_period
  ) %>%
  arrange(period, desc(percent_patients))

overall_yearly_pct

#=========== Summarize and create lineplot =====================================

label_map <- c("No AD, AP, or MS" = -1,
               "AD" = -3,
               "AP" = 2,
               "MS" = -1,
               "AD; AP" = -6,
               "AD; MS" = -2,
               "AP; MS" = 9,
               "AD; AP; MS" = 5)

# Inputs:
#   final_df: Final dataframe of drug combinations per patient and month
#   label_map: Named vector mapping drug combinations to label positions
#   box_padding: Numeric specifying how far away the labels should be
#   force_input: Numeric specifying amount of propulsion between labels
create_lineplot <- function(final_df, label_map, box_padding = 2,
                            force_input = 8) {
  # Number of patients
  n_patients <- length(unique(final_df$patid))
  
  # Summary table
  drug_summary_updated <- final_df %>%
    select(starts_with("month_")) %>%
    tbl_summary(
      by = NULL,
      statistic = list(all_categorical() ~ "{n} ({p}%)"),
      missing = "no"
    )
  
  # Convert wide to long
  long_summary <- final_df %>%
    pivot_longer(
      cols = -patid,
      names_to = "month",
      values_to = "drug_combination"
    )
  
  # Calculate prevalence of each drug combination
  # Denominator is all individuals in cohort (including those with no prescriptions)
  # contributing non-censored time. Individuals are dropped from the month 
  # after censoring
  summary_long <- long_summary %>%
    filter(!is.na(drug_combination)) %>%
    group_by(month, drug_combination) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(month) %>%
    mutate(percent = n / sum(n) * 100) %>%
    ungroup()
  
  # Add numeric month and relevel drug combination
  summary_long <- summary_long %>%
    mutate(
      month_num = as.integer(str_extract(month, "-?\\d+")),
      month_label = factor(month_num, levels = c(-12:-1, 1:12), 
                           labels = paste0("Month ", c(-12:-1, 1:12))),
      drug_combination = factor(drug_combination, 
                                levels = c("None", "AD", "AP", "MS", "AD; AP", 
                                           "AD; MS", "AP; MS", "AD; AP; MS"),
                                labels = c("No AD, AP, or MS", "AD", "AP", "MS", 
                                           "AD; AP", "AD; MS", "AP; MS", 
                                           "AD; AP; MS"))
    )
  
  # Default label map if no input specified
  if (is.null(label_map)) {
    label_map <- c("No AD, AP, or MS" = -1,
                   "AD" = -3,
                   "AP" = 2,
                   "MS" = -1,
                   "AD; AP" = -6,
                   "AD; MS" = -2,
                   "AP; MS" = 9,
                   "AD; AP; MS" = 5)
  }
  
  # Turn the named vector into case_match formulas: "AD" ~ -3, etc
  pairs <- imap(label_map, ~ new_formula(expr(!!.y), expr(!!.x)))
  
  # Make label positions staggered by drug combination
  summary_labels <- summary_long %>%
    group_by(drug_combination) %>%
    # no filtering: keep all months for lookup
    mutate(
      label_x = case_match(drug_combination, !!!pairs, 
                           # "No AD, AP, or MS" ~ -1,
                           # "AD" ~ -3,
                           # "AP" ~ 2,
                           # "MS" ~ -1,
                           # "AD; AP" ~ -6,
                           # "AD; MS" ~ -2,
                           # "AP; MS" ~ 9,
                           # "AD; AP; MS" ~ 5,
                           .default = 1
      ),
      label_y = percent[match(label_x, month_num)]
    ) %>%
    distinct(drug_combination, label_x, label_y, .keep_all = TRUE) %>%
    ungroup()
  
  # Plot (no title), colorblind-friendly, logical color combinations, different 
  # line types for none/single vs dual vs triple
  percentage_lineplot <- summary_long %>%
    ggplot(aes(x = month_num, y = percent, color = drug_combination, 
               linetype = drug_combination
    )) +
    annotate("rect",
             xmin = -1, xmax = 1,
             ymin = -Inf, ymax = Inf,
             fill = "grey90",
             alpha = 0.4) +
    geom_line(linewidth = 0.8, alpha = 0.85) +  # moderate transparency
    geom_point(size = 1.3, alpha = 0.85) +       # slightly smaller, semi-transparent points
    # Labels at last observed month for each combination
    geom_label_repel(
      data = summary_labels,
      aes(x = label_x, y = label_y, label = drug_combination),
      nudge_x = 0,
      nudge_y = 0,
      size = 3,
      segment.curvature = 0,
      segment.ncp = 0,
      segment.size = 0.2,
      show.legend = FALSE,
      min.segment.length = 0,
      box.padding = box_padding, # Increase: pushes labels further
      point.padding = 0.5, # Increase: pushes label further from point
      force = force_input, # More repulsion force among lables
      max.time = 8, # More iterations to find better spacing
      # segment.color = "gray70",
      max.overlaps = Inf
      # arrow = arrow(length = unit(0.02, "npc"))
    ) +
    scale_x_continuous(
      breaks = c(-12:12),
      labels = c(paste0("Month ", c(-12:-1)), "Diagnosis", paste0("Month ", c(1:12))),
      expand = expansion(add = c(0.2, 0.2)) # small spacing at plot edges
    ) +
    scale_y_continuous(
      limits = c(0, 60),
      breaks = seq(0, 60, 10),       # labels every 10%
      minor_breaks = seq(0, 60, 5),  # gridlines every 5%
      expand = expansion(mult = c(0, 0.05))
    ) +
    # Updated colorblind-friendly palette
    scale_color_manual(
      values = c(
        "AD" = "#E41A1C",          # bright red
        "AP" = "#377EB8",          # blue
        "MS" = "#FFD92F",          # yellow
        "AD; AP" = "#984EA3",      # purple/magenta
        "AD; MS" = "#FF7F00",      # orange
        "AP; MS" = "#4DAF4A",      # green
        "AD; AP; MS" = "#000000",  # black
        "No AD, AP, or MS" = "#999999"         # gray
      )
    ) +
    # Line types: solid (1 drug/None), dashed (2 drugs), dotted (3 drugs)
    scale_linetype_manual(
      values = c(
        "AD" = "solid",
        "AP" = "solid",
        "MS" = "solid",
        "AD; AP" = "dashed",
        "AD; MS" = "dashed",
        "AP; MS" = "dashed",
        "AD; AP; MS" = "dotted",
        "No AD, AP, or MS" = "solid"
      )
    ) +
    labs(
      x = "Month Relative to Diagnosis",
      y = "Percentage of Patients",
      color = "Drug Combination",
      linetype = "Drug Combination"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right", 
      panel.grid.minor.y = element_line(color = "gray95"), # new: light minor gridlines
      panel.grid.major.x = element_line(color = "gray90"), # vertical grid every month
      panel.grid.major.y = element_line(color = "gray90"), # horizontal grid every 10%
      panel.grid.minor.x = element_blank()
    ) + 
    annotate(
      "text",
      x = 11, y = 59,
      label = paste0("n = ", n_patients),
      hjust = 1, vjust = 1,              # anchors label at top-right *inside* panel
      size = 3.5
    )
  
  return(percentage_lineplot)
}

percentage_lineplot <- create_lineplot(final_df = final_df, 
                                       label_map = label_map)

percentage_lineplot

# # Save line plot
# ggsave(
#   plot = percentage_lineplot,
#   filename = paste0(wd, "Outputs/", "Percentage_lineplot", today(), ".png"),
#   dpi = 300, width = 12, height = 7, bg = "white")

# # Save as PDF
# ggsave(
#   plot = percentage_lineplot,
#   filename = paste0(wd, "Outputs/", "Percentage_lineplot", today(), ".pdf"),
#   dpi = 300, width = 12, height = 7, bg = "white")


#========== Separate line plots for pre- and post-covid ========================

# Separate out those diagnosed before 23 Mar 2020 (first UK stay-at-home order)
pre_covid <- bipolar_diagnosis_date %>% 
  filter(bipolar_diag_date < as.Date("2020-03-23")) 
final_df_pre_covid <- final_df %>%
  filter(patid %in% pre_covid$patid)

# Separate out those dignosed on/after COVID period (from 23 Mar 2020 and beyond)
from_covid <- bipolar_diagnosis_date %>% 
  filter(bipolar_diag_date >= as.Date("2020-03-23")) 
final_df_from_covid <- final_df %>%
  filter(patid %in% from_covid$patid)

# Create lineplots
percentage_lineplot_pre_covid <- create_lineplot(final_df = final_df_pre_covid,
                                                 label_map = label_map,
                                                 box_padding = 2,
                                                 force_input = 5)
percentage_lineplot_pre_covid

# # Save pre-covid line plot
# ggsave(
#   plot = percentage_lineplot_pre_covid,
#   filename = paste0(wd, "Outputs/", "Percentage_lineplot_pre_covid",
#                     today(), ".png"),
#   dpi = 300, width = 12, height = 7, bg = "white")

percentage_lineplot_from_covid <- create_lineplot(final_df = final_df_from_covid,
                                                  label_map = label_map,
                                                  box_padding = 2.6,
                                                  force_input = 2)
percentage_lineplot_from_covid

# # Save from-covid line plot
# ggsave(
#   plot = percentage_lineplot_from_covid,
#   filename = paste0(wd, "Outputs/", "Percentage_lineplot_from_covid",
#                     today(), ".png"),
#   dpi = 300, width = 12, height = 7, bg = "white")

# Combine the pre- and from- covid plots and save
combined <- ggarrange(
  percentage_lineplot_pre_covid + 
    ggtitle("Before COVID-19") +
    theme(plot.title = element_text(hjust = 0.5)), 
  percentage_lineplot_from_covid + 
    ggtitle("During and After COVID-19") +
    theme(plot.title = element_text(hjust = 0.5)), nrow = 2)
combined
# # Save combined plot
# ggsave(
#   plot = combined,
#   filename = paste0(wd, "Outputs/", "Percentage_lineplot_pre_from_covid_combined",
#                     today(), ".png"),
#   dpi = 300, width = 11, height = 9, bg = "white")

# # Save as PDF
# ggsave(
#   plot = combined,
#   filename = paste0(wd, "Outputs/", "Percentage_lineplot_pre_from_covid_combined",
#                     today(), ".pdf"),
#   dpi = 300, width = 11, height = 9, bg = "white")


#========== Monotherapy defined as single-agent ================================

## Single-agent summaries ##
# For each patient and each month since diagnosis, combine all unique drug 
# ingredients prescribed into a single semicolon-separated list
monthly_ingredients <- prescriptions_with_month %>%
  mutate(
    drug_ingredient = paste0(drug_type, "_", ingredient)
  ) %>%
  group_by(patid, month_since_diag) %>%
  summarise(
    ingredient_combined = paste(unique(drug_ingredient), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(patid, month_since_diag) 

# Pivot to wide format — one row per patid
wide_monthly_ingredients <- monthly_ingredients %>%
  pivot_wider(
    names_from = month_since_diag,
    values_from = ingredient_combined,
    names_prefix = "month_"
  )

# Add missing month columns (month_-12 to month_12) and replace NAs with "None"
all_months <- paste0("month_", c(-12:-1, 1:12))
wide_monthly_ingredients <- wide_monthly_ingredients %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(all_of(setdiff(all_months, names(.))), ~NA_character_)) %>%
  mutate(across(all_of(all_months), ~replace_na(.x, "None")))

# Incorporate prescription data into full cohort
final_df_monthly_ingredients <- unique_patids %>%
  left_join(wide_monthly_ingredients, by = "patid") %>%
  mutate(across(starts_with("month_"), ~replace_na(.x, "None")))

# Reorder columns
final_df_monthly_ingredients <- final_df_monthly_ingredients %>%
  select(patid, all_of(all_months))

# Incorporate censoring information into final dataset
# All months after censoring are set to NA
# Note: "None" means no prescription recorded, whereas NA means censored
final_df_ingredients <- final_df_monthly_ingredients %>%
  left_join(
    censor_dates %>% select(patid, month_since_diag, within_12_months),
    by = "patid"
  ) %>%
  mutate(across(all_of(paste0("month_", 1:12)), ~as.character(.x))) %>%
  rowwise() %>%
  mutate(across(
    paste0("month_", 1:12),
    ~if (!is.na(month_since_diag) && within_12_months &&
         as.integer(str_extract(cur_column(), "\\d+")) >= month_since_diag) {
      NA_character_
    } else {
      .x
    }
  )) %>%
  ungroup() %>%
  select(-month_since_diag, -within_12_months)

# Overall % of patients per drug combo pre- and post-diagnosis 
# Long format with month number and period (pre vs post)
long_year_ingredients <- final_df_ingredients %>%
  pivot_longer(
    cols = starts_with("month_"),
    names_to = "month",
    values_to = "ingredient_combination"
  ) %>%
  mutate(
    month_num = as.integer(str_extract(month, "-?\\d+")),
    period = case_when(
      month_num < 0 ~ "Year pre-diagnosis",
      month_num > 0 ~ "Year post-diagnosis",
      TRUE ~ NA_character_   # should drop the 0 month if it ever existed
    )
  ) %>%
  # Drop censored months (NA) and months outside the ±12 window (already handled)
  filter(!is.na(period), !is.na(ingredient_combination))


# Denominator: number of patients with any observed time in each period
denom_period_ingredients <- long_year_ingredients %>%
  distinct(patid, period) %>%
  count(period, name = "n_patients_period")

# For each period, count patients who ever had each drug combination at least once
overall_yearly_pct_ingredients <- long_year_ingredients %>%
  distinct(patid, period, ingredient_combination) %>%   # "ever on" that combo in that year
  count(period, ingredient_combination, name = "n_patients") %>%
  left_join(denom_period_ingredients, by = "period") %>%
  mutate(
    percent_patients = 100 * n_patients / n_patients_period
  ) %>%
  arrange(period, desc(percent_patients))

overall_yearly_pct_ingredients

# Number of patients
n_patients <- length(unique(final_df_ingredients$patid))

# Count the number of unique ingredients within each drug class, for each 
# patient and month, for year pre- and post-diagnosis
overall_yearly_ingredients <- overall_yearly_pct_ingredients %>%
  mutate(AD_times = str_count(ingredient_combination, "AD"),
         AP_times = str_count(ingredient_combination, "AP"),
         MS_times = str_count(ingredient_combination, "MS"),
         drug_combo = paste0("AD_", AD_times, "_AP_", AP_times, "_MS_", MS_times))

# Count the number of times that multiple drugs within the same drug class
# were prescribed, for year pre- and post-diagnosis
overall_yearly_combo <- overall_yearly_ingredients %>%
  group_by(period, drug_combo) %>%
  summarise(sum_n_patients = sum(n_patients), 
            n_patients_period = first(n_patients_period),
            ingredient_combos = paste(unique(ingredient_combination), 
                                      collapse = " OR ")) %>%
  mutate(sum_percent_patients = 100 * sum_n_patients / n_patients_period)

# Filter to pre-diagnosis
overall_yearly_combo_pre <- overall_yearly_combo %>%
  filter(period == "Year pre-diagnosis") %>%
  arrange(desc(sum_percent_patients))

# Filter to post-diagnosis
overall_yearly_combo_post <- overall_yearly_combo %>%
  filter(period == "Year post-diagnosis") %>%
  arrange(desc(sum_percent_patients))

## Filter to combinations with more than 1 drug within the same drug class
overall_yearly_combo_therapy_pre <- overall_yearly_combo_pre %>%
  filter(str_detect(drug_combo, "[23456789]"))
overall_yearly_combo_therapy_post <- overall_yearly_combo_post %>%
  filter(str_detect(drug_combo, "[23456789]"))

# How many patients had multiple drugs within the same class prescribed for 
# at least one month in the year pre- or post-diagnosis: 23.7% pre, 25.5% post
sum(overall_yearly_combo_therapy_pre$sum_percent_patients)
sum(overall_yearly_combo_therapy_post$sum_percent_patients)


### Create table of multiple drugs within the same class among those with 
### monotherapy

## AD class monotherapy
# Filter to those with AD monotherapy for each month
AD_monotherapy <- long_year_ingredients %>%
  filter(str_count(ingredient_combination, "AD") >= 1 & 
           str_count(ingredient_combination, "AP") == 0 & 
           str_count(ingredient_combination, "MS") == 0)
# Create variable for multiple drugs within the same class
AD_monotherapy <- AD_monotherapy %>%
  mutate(count_multiple = str_count(ingredient_combination, "AD")) %>%
  # combine 3+ into one category
  mutate(count_multiple = ifelse(count_multiple >= 3, "3+", count_multiple))

# Calculate prevalence of multiple ingredients within AD monotherapy
# Denominator is all individuals on AD class monotherapy in each month
# contributing non-censored time. Individuals are dropped from the month 
# after censoring
summary_AD_monotherapy <- AD_monotherapy %>%
  group_by(month, count_multiple) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(month) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()
# Summary table
summary_AD_monotherapy_table <- summary_AD_monotherapy %>%
  select(starts_with("month_")) %>%
  tbl_summary(
    by = NULL,
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    missing = "no"
  )

## AP class monotherapy
# Filter to those with AP monotherapy for each month
AP_monotherapy <- long_year_ingredients %>%
  filter(str_count(ingredient_combination, "AP") >= 1 & 
           str_count(ingredient_combination, "AD") == 0 & 
           str_count(ingredient_combination, "MS") == 0)
# Create variable for multiple drugs within the same class
AP_monotherapy <- AP_monotherapy %>%
  mutate(count_multiple = str_count(ingredient_combination, "AP")) %>%
  # combine 3+ into one category
  mutate(count_multiple = ifelse(count_multiple >= 3, "3+", count_multiple))
# Calculate prevalence of multiple ingredients within AP monotherapy
summary_AP_monotherapy <- AP_monotherapy %>%
  group_by(month, count_multiple) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(month) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()


## MS class monotherapy
# Filter to those with MS monotherapy for each month
MS_monotherapy <- long_year_ingredients %>%
  filter(str_count(ingredient_combination, "MS") >= 1 & 
           str_count(ingredient_combination, "AP") == 0 & 
           str_count(ingredient_combination, "AD") == 0)
# Create variable for multiple drugs within the same class
MS_monotherapy <- MS_monotherapy %>%
  mutate(count_multiple = str_count(ingredient_combination, "MS")) %>%
  # combine 3+ into one category
  mutate(count_multiple = ifelse(count_multiple >= 3, "3+", count_multiple))
# Calculate prevalence of multiple ingredients within MS monotherapy
summary_MS_monotherapy <- MS_monotherapy %>%
  group_by(month, count_multiple) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(month) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()


### Helper function to create table examining within-class multiple drugs
# Format the mutliple drugs within the same class monotherapy tables
# df_sub: input drug table to format
# group_name: name of monotherapy group
# month_levels: vector of month levels
# digits: number of digits to round to
make_med_table <- function(df_sub, group_name, month_levels, digits = 1) {
  # Main 3-row table: 1, 2, 3+
  main_part <- df_sub %>%
    mutate(
      count_multiple = factor(count_multiple, levels = c("1", "2", "3+")),
      month_label = paste("Month", sub("month_", "", month)),
      cell = paste0(n, "\n(", format(round(percent, digits), digits), "%)")
    ) %>%
    select(count_multiple, month_label, cell) %>%
    distinct() %>%
    complete(
      count_multiple,
      month_label = month_levels,
      fill = list(cell = "")
    ) %>%
    pivot_wider(
      id_cols = count_multiple,
      names_from = month_label,
      values_from = cell
    ) %>%
    arrange(count_multiple) %>%
    mutate(
      `Number of Medications` = as.character(count_multiple),
      group = group_name
    ) %>%
    select(group, `Number of Medications`, all_of(month_levels))
  
  # Total row: sum n across medication categories for each month
  total_part <- df_sub %>%
    mutate(
      month_label = paste("Month", sub("month_", "", month))
    ) %>%
    group_by(month_label) %>%
    summarise(total_n = sum(n, na.rm = TRUE), .groups = "drop") %>%
    complete(
      month_label = month_levels,
      fill = list(total_n = 0)
    ) %>%
    mutate(
      total_n = as.character(total_n),
      group = group_name,
      `Number of Medications` = "Total"
    ) %>%
    select(group, `Number of Medications`, month_label, total_n) %>%
    pivot_wider(
      id_cols = c(group, `Number of Medications`),
      names_from = month_label,
      values_from = total_n
    ) %>%
    select(group, `Number of Medications`, all_of(month_levels))
  
  bind_rows(main_part, total_part)
}

### Create pre-diagnosis tables for each drug class, then combine

month_levels_pre <- paste("Month", -12:-1)
tab_AD_pre <- make_med_table(summary_AD_monotherapy, 
                             "Antidepressant Class Monotherapy", 
                             month_levels = month_levels_pre)
tab_AP_pre <- make_med_table(summary_AP_monotherapy, 
                             "Antipsychotic Class Monotherapy", 
                             month_levels = month_levels_pre)
tab_MS_pre <- make_med_table(summary_MS_monotherapy, 
                             "Mood Stabiliser Class Monotherapy", 
                             month_levels = month_levels_pre)
combined_tab_pre <- bind_rows(tab_AD_pre, tab_AP_pre, tab_MS_pre) %>%
  as_grouped_data(
    groups = "group",
    columns = c("Number of Medications", month_levels_pre)
  )

# Display pre-diagnosis table
# gt_table_pre <- combined_tab_pre %>%
#   gt(groupname_col = "group") %>%
#   cols_hide(columns = group) %>%
#   cols_align(
#     align = "left",
#     columns = `Number of Medications`
#   ) %>%
#   cols_align(
#     align = "center",
#     columns = -`Number of Medications`
#   ) %>%
#   tab_style(
#     style = cell_text(weight = "bold"),
#     locations = cells_column_labels(everything())
#   )
# gt_table_pre

ft_pre <- combined_tab_pre %>%
  flextable() %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%          # smaller text
  padding(padding = 2, part = "all") %>%        # tighter cells
  set_table_properties(
    layout = "autofit",
    width = 1
  )


### Create post-diagnosis tables for each drug class, then combine

month_levels_post <- paste("Month", 1:12)
tab_AD_post <- make_med_table(summary_AD_monotherapy, 
                             "Antidepressant Class Monotherapy", 
                             month_levels = month_levels_post)
tab_AP_post <- make_med_table(summary_AP_monotherapy, 
                             "Antipsychotic Class Monotherapy", 
                             month_levels = month_levels_post)
tab_MS_post <- make_med_table(summary_MS_monotherapy, 
                             "Mood Stabiliser Class Monotherapy", 
                             month_levels = month_levels_post)
combined_tab_post <- bind_rows(tab_AD_post, tab_AP_post, tab_MS_post) %>%
  as_grouped_data(
    groups = "group",
    columns = c("Number of Medications", month_levels_post)
  )

# Display post-diagnosis table
ft_post <- combined_tab_post %>%
  flextable() %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%          # smaller text
  padding(padding = 2, part = "all") %>%        # tighter cells
  set_table_properties(
    layout = "autofit",
    width = 1
  )

landscape_section <- prop_section(
  page_size = page_size(
    orient = "landscape",
    width = 11.7,
    height = 8.3
  ), 
  type = "continuous",
  page_margins = page_mar()
)
## Save pre-diagnosis and post-diagnosis tables to word document
doc <- read_docx() %>%
  body_end_section_continuous() %>%
  body_set_default_section(landscape_section) %>%
  body_add_flextable(ft_pre) %>%
  body_add_par("") %>%
  body_add_flextable(ft_post)
print(doc, target = paste0(wd, "Outputs/", "multiple_meds_within_class", 
                           today(), ".docx"))


#========== Miscellaneous ======================================================

# Rate of BP diagnoses over time
bipolar_cohort %>% 
  count(bipolar_year) %>%
  ggplot(aes(x = bipolar_year, y = n)) +
  geom_col(fill = "#377EB8", alpha = 0.8) +
  geom_text(aes(label = n), vjust = -0.4, size = 3) +
  scale_x_continuous(
    breaks = seq(min(bipolar_cohort$bipolar_year, na.rm = TRUE),
                 max(bipolar_cohort$bipolar_year, na.rm = TRUE),
                 by = 1)  # one tick per year
  ) +
  labs(
    title = "Number of People Diagnosed with Bipolar Disorder by Year",
    x = "Year of Diagnosis",
    y = "Number of Diagnoses"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )




