# -----------------------
# BIPOLAR DISORDER PRESCRIBING TRENDS 2000-2022
# -----------------------
# CPRD 2023 data

# Lithium and AD monotherapy regression analyses (complete case and multiple imputation) ####
# Last updated: 23/12/2025
# 
# Data:
#   - Data/bipolar_cohort_ads.rdata: Antidepressant prescriptions
#   - Data/bipolar_cohort_moodstabs.rdata: Mood stabiliser prescriptions
#   - Data/bipolar_cohort_aps.rdata: Antipsychotic prescriptions
#   - Data/bipolar_cohort.rdata: Patient information 
# 
# Outputs:
#   - Outputs/lithium_regression_combined: Regression output for lithium models  
#   - Outputs/ad_regression_combined: Regression output for antidepressant models
#   - Outputs/lithiumAD_forest_plot: Forest plot of lithium and antidepressant regression output
# -----------------------


# Clear memory
rm(list = ls())

# Packages
library(lubridate)
library(gtsummary)
library(tidyverse)
library(purrr)
library(gt)
library(tidylog)
library(data.table)
library(mice)
library(finalfit)
library(dplyr)
library(forestplot)
library(tidyr)
library(grid)
library(janitor)
library(ggplot2)


# Get working directory
wd <- getwd()

# Set output directory
output_dir <- "Outputs/"

# Load data
load("Data/bipolar_cohort_ads.rdata")
load("Data/bipolar_cohort_moodstabs.rdata")
load("Data/bipolar_cohort_aps.rdata")
load("Data/bipolar_cohort.rdata")

# Data management

# Censor at earliest of death, lcd, or regendate -- include gender
censor_dates <- bipolar_cohort %>%
  select(patid, gender, deathdate, regenddate, lcd, first_bipolar_date) %>%
  mutate(deathdate = as.Date(deathdate),
         regenddate = as.Date(regenddate),
         lcd = as.Date(lcd),
         earliest_exit = pmin(deathdate, regenddate, lcd, na.rm = TRUE),
         diff_since_diag = earliest_exit - first_bipolar_date,
         month_since_diag = ifelse(diff_since_diag >= 0,  # (should all be >= 0)
                                   ceiling(diff_since_diag / 30.42),  # after dianosis
                                   floor(diff_since_diag / 30.42)),   # before diagnosis
         within_12_months = !is.na(month_since_diag) & abs(month_since_diag) <= 12) %>%
  select(patid, gender, earliest_exit, diff_since_diag, month_since_diag, within_12_months)

# Prescription data

# Combine prescriptions from all drug classes
ads <- bipolar_cohort_ads %>%
  select(patid, issuedate, ingredient = Antidepressant) %>%
  mutate(drug_type = "AD")

aps <- bipolar_cohort_aps %>%
  select(patid, issuedate, ingredient = AP) %>%
  mutate(drug_type = "AP")

moodstabs <- bipolar_cohort_moodstabs %>%
  select(patid, issuedate, ingredient = MoodStabiliser) %>%
  mutate(drug_type = "MS")

# Combine all prescriptions & remove duplicates
all_prescriptions <- bind_rows(ads, aps, moodstabs) %>%
  distinct() %>%
  inner_join(select(bipolar_cohort, patid, first_bipolar_date), by = "patid") %>%
  mutate(diff_days = as.numeric(issuedate - first_bipolar_date)) %>%
  filter(diff_days >= 0 & diff_days <= 365) %>%
  mutate(cutoff_date = first_bipolar_date + 365) %>%
  left_join(censor_dates %>% select(patid, gender, earliest_exit), by = "patid")

rm(ads, aps, moodstabs, bipolar_cohort_ads, bipolar_cohort_aps, bipolar_cohort_moodstabs)

# n_distinct(all_prescriptions$patid) # 37,208

#=============== Get drug combinations in treatment intervals

# Parameters
overlap_days <- 30  # drugs within 30 days are considered combination
interrupt_cutoff_months <- 3  # >= 3 consecutive months

get_treatment_intervals <- function(prescriptions_prep, overlap_days, 
                                    interrupt_cutoff_months) {
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
      ingredient_combination = {
        ing <- sort(unique(na.omit(ingredient)))
        if (length(ing) == 0) {
          NA_character_
        } else {
          paste(ing, collapse = "; ")
        }
      },
      int_start = first(int_start),
      int_end = first(int_end),
      first_rx_date = first(first_rx_date),
      earliest_exit = first(earliest_exit),
      gender = first(gender),
      .groups = "drop") %>%
    ungroup()

  # If censored day lies within a time interval, all intervals after are censored. 
  # Filter out those with censored date before int_start
  prescriptions_combs_censor <- prescriptions_combs %>%
    mutate(drug_combination = ifelse(earliest_exit < int_start, "Censored", 
                                     drug_combination),
           ingredient_combination = ifelse(earliest_exit < int_start, "Censored",
                                           ingredient_combination)) %>%
    filter(earliest_exit > int_start) %>%
    group_by(patid) %>%
    # Get maximum interval for each patient
    mutate(max_int_idx = max(int_idx, na.rm = TRUE)) %>%
    ungroup()
  
  # Fill in intervals with missing prescription
  prescriptions_combs_fill <- prescriptions_combs_censor %>%
    group_by(patid) %>%
    complete(int_idx = full_seq(1:unique(max_int_idx), 1)) %>%
    mutate(drug_combination = replace_na(drug_combination, "None"),
           ingredient_combination = replace_na(ingredient_combination, "None")) %>%
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
        is_none & (run_length >= interrupt_cutoff_months) ~ "Int/Disc",
        TRUE ~ drug_combination),
      ingredient_combination = case_when(
        is_none & (run_length >= interrupt_cutoff_months) ~ "Int/Disc",
        TRUE ~ ingredient_combination))
  
  return(prescriptions_combs_final)}

# Run function
prescriptions_combs_final <- get_treatment_intervals(
  prescriptions_prep = all_prescriptions, overlap_days = overlap_days, 
  interrupt_cutoff_months = interrupt_cutoff_months)

# Create binary indicators at patient level
patient_lithium_ad <- prescriptions_combs_final %>%
  group_by(patid) %>%
  summarise(lithium_prescribed = as.integer(any(grepl("Lithium", ingredient_combination, ignore.case = TRUE), na.rm = TRUE)),
            ad_only_prescribed = as.integer(any(drug_combination == "AD", na.rm = TRUE)),
            .groups = "drop")

 n_distinct(patient_lithium_ad$patid) # 37,067 - a small number (141) lost due to incomplete treatment intervals
 n_distinct(all_prescriptions$patid) # 37,208
 
incomplete_intervals <- all_prescriptions %>%
   select(patid) %>%
   distinct() %>%
   anti_join(patient_lithium_ad, by = "patid") %>%
   mutate(incomplete = 1)
 
 rm(all_prescriptions, prescriptions_combs_final)

#========= Combine prescription and cohort data 
 
model_data <- bipolar_cohort %>%
   select(patid, gender, ethnicity_cat_cprdhes, pat_2019imd_quintile, bipolar_year, age_at_first_bipolar) %>%
   mutate(
     ethnicity_cat_cprdhes = str_to_title(ethnicity_cat_cprdhes),
     ethnicity_cat_cprdhes = case_when(
       ethnicity_cat_cprdhes == "Mixed" ~ "Mixed/Other",
       ethnicity_cat_cprdhes == "Other" ~ "Mixed/Other",
       TRUE ~ ethnicity_cat_cprdhes),
     ethnicity_cat_cprdhes = factor(ethnicity_cat_cprdhes,
                                    levels = c("White", "Black", "Asian", "Mixed/Other"))) %>%
   full_join(patient_lithium_ad, by = "patid") %>%
   mutate( # recoding missing to 0 (people that did not receive prescriptions within 1y)
     lithium_prescribed = if_else(is.na(lithium_prescribed), 0L, lithium_prescribed),
     ad_only_prescribed = if_else(is.na(ad_only_prescribed), 0L, ad_only_prescribed),
     time_period = case_when(bipolar_year < 2006 ~ "2000–2005", # Derive time period variable
                                    bipolar_year >= 2006 & bipolar_year <= 2010 ~ "2006–2010",
                                    bipolar_year >= 2011 & bipolar_year <= 2015 ~ "2011–2015",
                                    bipolar_year >= 2016 & bipolar_year <= 2019 ~ "2016–2019",
                                    bipolar_year >= 2020 ~ "2020–2022")) %>%
  left_join(incomplete_intervals, by = "patid") %>%
  mutate(across(
    c(lithium_prescribed, ad_only_prescribed),
    ~ if_else(incomplete %in% 1, NA_integer_, .))) %>% # set those with missing intervals to NA
  select(-incomplete)

rm(incomplete_intervals)
     
n_distinct(model_data$patid) # 40,965

# =============================================================================
# Testing linearity of diagnosis year x outcome
# =============================================================================

plot_year_linearity <- function(data, outcome_var, outcome_label, title_prefix = NULL) {
  
  # Compute yearly proportions
  plot_df <- data %>%
    group_by(bipolar_year) %>%
    summarise(
      prop = mean(.data[[outcome_var]], na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  # --- Fit linear vs quadratic models for p-value ---
  m_lin  <- lm(prop ~ bipolar_year, data = plot_df)
  m_quad <- lm(prop ~ poly(bipolar_year, 2, raw = TRUE), data = plot_df)
  anova_p <- anova(m_lin, m_quad)$`Pr(>F)`[2]
  p_text_val <- ifelse(anova_p < 0.001, "<0.001", round(anova_p, 2))
  p_text <- paste0("ANOVA p (Linear vs Linear + Quadratic): ", p_text_val)
  
  # Define smoother specs
  smoothers <- tribble(
    ~label,               ~method, ~formula,                   ~colour,    ~fill,                       ~linetype,
    "Linear",             "lm",    NULL,                       "#D55E00", alpha("#D55E00", 0.25), "dashed",
    "Linear + Quadratic", "lm",    y ~ poly(x, 2, raw = TRUE), "#009E73", alpha("#009E73", 0.25), "dotdash"
  )
  
  # Base plot
  p <- ggplot(plot_df, aes(x = bipolar_year, y = prop)) +
    geom_line(aes(group = 1), colour = "black", linewidth = 1, alpha = 0.7) +
    geom_point(aes(size = n), alpha = 0.7)
  
  # Add smoothers
  walk(seq_len(nrow(smoothers)), \(i) {
    p <<- p + geom_smooth(
      aes(
        colour   = smoothers$label[i],
        fill     = smoothers$label[i],
        linetype = smoothers$label[i]
      ),
      method    = smoothers$method[i],
      formula   = smoothers$formula[[i]],
      se        = TRUE,
      linewidth = 0.9
    )
  })
  
  # Scales, labels, theme
  p <- p +
    scale_colour_manual(
      "Model type",
      values = deframe(smoothers[, c("label", "colour")])
    ) +
    scale_fill_manual(
      "Model type",
      values = deframe(smoothers[, c("label", "fill")])
    ) +
    scale_linetype_manual(
      "Model type",
      values = deframe(smoothers[, c("label", "linetype")])
    ) +
    scale_y_continuous(
      limits = c(0, 0.4)
    ) +
    labs(
      x = "Year of diagnosis",
      y = paste0("Proportion prescribed ", outcome_label)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
  
  # Conditionally add title
  if (!is.null(title_prefix)) {
    p <- p + labs(title = paste0(title_prefix, " by year of bipolar disorder diagnosis"))
  }
  
  p
}

# Generate plots (no titles)
plots <- list(
  lithium = plot_year_linearity(
    model_data,
    "lithium_prescribed",
    "lithium"
  ),
  ad_mono = plot_year_linearity(
    model_data,
    "ad_only_prescribed",
    "Antidepressant monotherapy"
  )
)

# Save plots
walk(names(plots), \(name) {
  ggsave(
    filename = file.path(
      output_dir,
      paste0("lineality_plot_", name, "_year_", today(), ".pdf")
    ),
    plot   = plots[[name]],
    width  = 8,
    height = 5,
    units  = "in"
  )
})

rm(plots)

# =============================================================================
# SETUP: Define variables and labels
# =============================================================================

all_vars <- c("gender", "pat_2019imd_quintile", "ethnicity_cat_cprdhes",
              "age_at_first_bipolar", "time_period")

var_labels <- list(
  gender = "Sex",
  pat_2019imd_quintile = "IMD quintile",
  ethnicity_cat_cprdhes = "Ethnicity",
  age_at_first_bipolar = "Age at first diagnosis",
  time_period = "Diagnosis time-period")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Define consistent color palette for characteristics
CHARACTERISTIC_COLORS <- c(
  "Sex" = "#2C7A7B",
  "IMD quintile" = "#C05621",
  "Ethnicity" = "#5A67D8",
  "Age at first diagnosis" = "#B83280",
  "Diagnosis time-period" = "#4C7C3E"
)

# =============================================================================
# MODEL FITTING FUNCTIONS
# =============================================================================

fit_models <- function(outcome, var, data = NULL, imp_data = NULL) {
  # Unified function to fit both CC and MI models
  # Args:
  #   outcome: outcome variable name
  #   var: predictor variable name
  #   data: complete case data (for CC models)
  #   imp_data: imputed data (for MI models)
  # Returns:
  #   List with unadjusted and adjusted model objects
  
  var <- as.character(var)
  
  if (!is.null(data)) {
    # Complete case models
    unadj <- glm(as.formula(paste(outcome, "~", var)),
                 data = data, family = "binomial")
    adj <- glm(as.formula(paste(outcome, "~", var, "+ bipolar_year + I(bipolar_year^2)")),
               data = data, family = "binomial")
  } else if (!is.null(imp_data)) {
    # Multiple imputation models
    unadj <- with(imp_data,
                  glm(as.formula(paste(outcome, "~", var)), family = "binomial"))
    adj <- with(imp_data,
                glm(as.formula(paste(outcome, "~", var, "+ bipolar_year + I(bipolar_year^2)")),
                    family = "binomial"))
  } else {
    stop("Must provide either data or imp_data")
  }
  
  return(list(unadj = unadj, adj = adj))
}

# =============================================================================
# COMBINED TABLE GENERATION
# =============================================================================

generate_combined_table <- function(outcome, cc_data, imp_data, vars_to_model) {
  cat("\nGenerating combined table for:", outcome, "\n")
  
  # Prepare mutually adjusted formula for MI (exclude time_period)
  mutual_vars <- setdiff(vars_to_model, "time_period")
  mutual_formula <- paste(outcome, "~", paste(c(mutual_vars, "bipolar_year", "I(bipolar_year^2)"), collapse = " + "))
  cat("Mutually adjusted formula for MI (excluding time_period):\n", mutual_formula, "\n\n")
  mi_mutual_model <- with(imp_data, glm(as.formula(mutual_formula), family = "binomial"))
  
  # Helper: create OR (95% CI) table
  make_combined_tbl <- function(model, include_var) {
    as.data.frame(
      tbl_regression(model, exponentiate = TRUE, include = include_var,
                     label = setNames(var_labels[include_var], include_var))
    ) %>%
      clean_names() %>%
      separate(x95_percent_ci, into = c("ci_low", "ci_high"), sep = ", ") %>%
      mutate(across(c(or, ci_low, ci_high), as.numeric),
             result = ifelse(is.na(or), NA,
                             sprintf("%.2f (%.2f, %.2f)", or, ci_low, ci_high))) %>%
      select(group = characteristic, result)
  }
  
  all_results <- lapply(vars_to_model, function(var) {
    cat("Processing variable:", var, "-", var_labels[[var]], "\n")
    var_char <- as.character(var)
    
    tryCatch({
      # Fit models
      cc_models <- fit_models(outcome, var_char, data = cc_data)
      mi_models <- fit_models(outcome, var_char, imp_data = imp_data)
      
      # Generate tables
      tbls <- list(
        cc_unadj = make_combined_tbl(cc_models$unadj, var_char),
        cc_adj = make_combined_tbl(cc_models$adj, var_char),
        mi_unadj = make_combined_tbl(mi_models$unadj, var_char),
        mi_adj = make_combined_tbl(mi_models$adj, var_char),
        mi_mutual = if (var == "time_period") {
          data.frame(group = unique(make_combined_tbl(mi_models$adj, var_char)$group),
                     result = NA_character_)
        } else {
          make_combined_tbl(mi_mutual_model, var_char)
        }
      )
      
      # Merge all tables
      Reduce(function(x, y) left_join(x, y, by = "group"), tbls) %>%
        setNames(c("group", "result_cc_unadj", "result_cc_adj", 
                   "result_mi_unadj", "result_mi_adj", "result_mi_mutual"))
      
    }, error = function(e) {
      cat("Error processing", var, ":", e$message, "\n")
      NULL
    })
  })
  
  cat("Stacking all variables...\n")
  
  # Calculate missingness for footnote
  # --- Missingness calculation ---
  # Include outcome in the calculation
  all_vars_check <- c(outcome, vars_to_model)
  
  missing_pcts <- sapply(all_vars_check, function(v) mean(is.na(cc_data[[v]])) * 100)
  
  # Separate fully observed vs variables with missing data
  fully_observed <- names(missing_pcts)[missing_pcts == 0]
  has_missing <- names(missing_pcts)[missing_pcts > 0]
  
  # Build text for fully observed
  fully_obs_text <- paste(sapply(fully_observed, function(v) {
    if (v == outcome) return("the outcome")
    var_labels[[v]]
  }), collapse = ", ")
  
  # Build text for variables with missing data
  missing_text <- paste(sapply(has_missing, function(v) {
    pct <- missing_pcts[v]
    pct_str <- if (pct > 0 & pct < 0.05) "<0.1%" else sprintf("%.1f%%", pct)
    
    if (v == outcome) {
      sprintf("the outcome (%s)", pct_str)
    } else {
      sprintf("%s (%s)", var_labels[[v]], pct_str)
    }
  }), collapse = ", ")
  
  # Combine into final note
  missingness_note <- sprintf(
    "Fully observed variables: %s; variables with missing data: %s",
    fully_obs_text, missing_text
  )
  
  # Print for verification
  cat("Missingness note:", missingness_note, "\n")
  
  # Combine and format
  combined_df <- do.call(rbind, all_results) %>%
    rename(`CC Unadjusted` = result_cc_unadj,
           `CC Time-adjusted` = result_cc_adj,
           `MI Unadjusted` = result_mi_unadj,
           `MI Time-adjusted` = result_mi_adj,
           `MI Multivariable Model` = result_mi_mutual) %>%
    mutate(across(2:6, ~ ifelse(trimws(.) == "NA (NA, NA)", NA, .)))
  
  # Identify heading rows and blank estimates (except Age)
  heading_rows <- combined_df$group %in% unlist(var_labels)
  combined_df <- combined_df %>%
    mutate(across(2:6, ~ ifelse(heading_rows & !grepl("Age", group, ignore.case = TRUE), "", .)))
  
  # Create gt table
  combined_df %>%
    gt(rowname_col = "group") %>%
    tab_stubhead(label = "Characteristic") %>%
    sub_missing(missing_text = "-") %>%
    cols_align(align = "center", columns = 2:6) %>%
    tab_style(style = cell_text(weight = "bold"),
              locations = list(cells_stubhead(), 
                               cells_stub(rows = heading_rows),
                               cells_column_labels())) %>%
    tab_style(style = cell_text(indent = px(20)),
              locations = cells_stub(rows = !heading_rows)) %>%
    tab_footnote("CC = Complete case; univariate logistic regression (unadjusted)",
                 cells_column_labels("CC Unadjusted")) %>%
    tab_footnote("CC = Complete case; logistic regression adjusted for diagnosis year and diagnosis year squared",
                 cells_column_labels("CC Time-adjusted")) %>%
    tab_footnote("MI = Multiple imputation; univariate logistic regression (unadjusted)",
                 cells_column_labels("MI Unadjusted")) %>%
    tab_footnote("MI = Multiple imputation; logistic regression adjusted for diagnosis year and diagnosis year squared",
                 cells_column_labels("MI Time-adjusted")) %>%
    tab_footnote("MI = Multiple imputation; multivariable logistic regression including all covariates (excluding time period) plus diagnosis year and diagnosis year squared",
                 cells_column_labels("MI Multivariable Model")) %>%
    tab_source_note("CC = Complete case; MI = Multiple imputation; IMD = Index of Multiple Deprivation") %>%
    tab_source_note("All estimates are odds ratios with 95% confidence intervals") %>%
    tab_source_note(missingness_note)}

# =============================================================================
# FOREST PLOT FUNCTIONS
# =============================================================================

REFERENCE_CHARS <- c( "Sex", "IMD quintile", "Ethnicity", "Age at first diagnosis", "Diagnosis time-period" )

assign_colors <- function(df) { df %>% mutate( char_group = 
                                                 case_when( characteristic %in% REFERENCE_CHARS ~ characteristic, TRUE ~ NA_character_ ) ) %>% 
    fill(char_group, .direction = "down") %>% 
    mutate( color = CHARACTERISTIC_COLORS[char_group], color = ifelse(is_header, paste0(color, "B3"), color) ) %>% select(-char_group) }

generate_forest_plot <- function(outcome, model_data_imp, vars_to_model, var_labels,
                                 digits_plot = 4, plot_title = NULL, add_footnote = TRUE,
                                 plot_max = NULL) {
  
  # Header labels mapping
  HEADER_REF_LABELS <- c(
    "Sex" = "Sex (Reference: Female)",
    "IMD quintile" = "IMD quintile (Reference: 1 [Least deprived])",
    "Ethnicity" = "Ethnicity (Reference: White)",
    "Age at first diagnosis" = "Age at first diagnosis",
    "Diagnosis time-period" = "Diagnosis time-period (Reference: 2000–2005)"
  )
  
  # Generate MI table
  all_results <- list()
  for (var in vars_to_model) {
    label_list <- list()
    label_list[[var]] <- var_labels[[var]]
    
    mi_models <- fit_models(outcome, var, imp_data = model_data_imp)
    mi_tbl <- if (var == "time_period") mi_models$unadj else mi_models$adj
    
    tbl_df <- as.data.frame(
      tbl_regression(
        mi_tbl,
        exponentiate = TRUE,
        label = label_list,
        estimate_fun = purrr::partial(style_ratio, digits = digits_plot)
      )
    ) %>%
      clean_names() %>%
      separate(x95_percent_ci, into = c("conf_low", "conf_high"), sep = ", ") %>%
      mutate(
        across(c(or, conf_low, conf_high), as.numeric),
        is_header = characteristic %in% unname(var_labels),
        plot_label = ifelse(is_header, HEADER_REF_LABELS[characteristic], paste0("  ", characteristic)),
        # Create CI label for all rows that have an estimate, including continuous vars
        ci_label = ifelse(!is.na(or),
                          paste0(format(round(or, 2), nsmall = 2),
                                 " (", format(round(conf_low, 2), nsmall = 2),
                                 "-", format(round(conf_high, 2), nsmall = 2), ")"),
                          "")
      ) %>%
      filter(!grepl("bipolar_year", characteristic)) %>%
      select(plot_label, ci_label, estimate = or, conf.low = conf_low,
             conf.high = conf_high, is_header, characteristic)
    
    all_results[[var]] <- tbl_df
  }
  
  mi_table <- do.call(rbind, all_results) %>%
    filter(!(is_header == FALSE & is.na(estimate)))  # remove reference rows
  
  # Prepare plot data
  plot_data <- mi_table %>%
    assign_colors() %>%
    mutate(
      row_num = rev(seq_len(n())),
      text_face = ifelse(is_header, "bold", "plain"),
      text_size = ifelse(is_header, 3.5, 3.2)
    )
  
  # Data-driven x-axis limits
  x_min <- min(plot_data$conf.low, na.rm = TRUE)
  x_max <- max(plot_data$conf.high, na.rm = TRUE)
  plot_min <- max(0, x_min - 0.06)
  if (is.null(plot_max)) {
    plot_max <- x_max + 0.05
  }
  arrow_threshold <- plot_max - 0.05
  
  # Legend
  legend_data <- data.frame(
    characteristic = names(CHARACTERISTIC_COLORS),
    color = unname(CHARACTERISTIC_COLORS),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      y = seq(0.95, 0.75, length.out = n()),
      label = case_when(
        characteristic == "Sex" ~ "Sex",
        characteristic == "IMD quintile" ~ "Deprivation",
        characteristic == "Ethnicity" ~ "Ethnicity",
        characteristic == "Age at first diagnosis" ~ "Age at diagnosis",
        characteristic == "Diagnosis time-period" ~ "Time period*",
        TRUE ~ characteristic
      )
    )
  
  # Build plot
  forest_plot <- ggplot(plot_data, aes(y = row_num)) +
    geom_segment(data = filter(plot_data, !is.na(estimate) & conf.high <= arrow_threshold),
                 aes(x = conf.low, xend = conf.high, yend = row_num, color = color), linewidth = 1) +
    geom_segment(data = filter(plot_data, !is.na(estimate) & conf.high > arrow_threshold),
                 aes(x = conf.low, xend = arrow_threshold, yend = row_num, color = color),
                 linewidth = 1,
                 arrow = arrow(angle = 20, length = unit(0.15, "cm"), type = "closed")) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.6) +
    geom_point(data = filter(plot_data, !is.na(estimate)),
               aes(x = estimate, fill = color), color = "white", size = 5, shape = 22, stroke = 1.5) +
    annotate("text", x = plot_min, y = plot_data$row_num,
             label = plot_data$plot_label, hjust = 0,
             size = plot_data$text_size, fontface = plot_data$text_face) +
    annotate("text", x = plot_max, y = plot_data$row_num,
             label = plot_data$ci_label, hjust = 1, size = 3) +
    scale_x_continuous(name = "Odds Ratio (95% CI)", limits = c(plot_min, plot_max), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0.015, 0.015)) +
    scale_color_identity() +
    scale_fill_identity() +
    coord_cartesian(xlim = c(plot_min, plot_max), clip = "off") +
    theme_minimal() +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(color = "gray92", linewidth = 0.5),
          axis.line.x = element_line(color = "black", linewidth = 0.6),
          axis.ticks.x = element_line(color = "black", linewidth = 0.5),
          axis.ticks.length.x = unit(0.15, "cm"),
          plot.margin = unit(c(0.5, 5, 0.8, 5.5), "cm"),
          axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 10)),
          axis.text.x = element_text(size = 10, color = "black"),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 15)),
          plot.caption = element_text(size = 9, hjust = 0, margin = margin(t = 10))) +
    ggtitle(plot_title)
  
  if (add_footnote) {
    forest_plot <- forest_plot + 
      labs(caption = paste0(
        "*All estimates are odds ratios (ORs) adjusted for diagnosis year and year²,\n",
        "except for time period which is shown unadjusted due to collinearity with year of diagnosis.\n",
        "Missing values were replaced using multiple imputation."
      ))
  }  
  return(forest_plot)}

# =============================================================================
# MAIN ANALYSIS PIPELINE
# =============================================================================

# Prepare and impute data
imputationvariables <- model_data %>%
  mutate_if(is.character, as.factor)

predM <- finalfit::missing_predictorMatrix(
  imputationvariables,
  drop_from_imputed = c("patid"),
  drop_from_imputer = c("patid"))

model_data_imp <- futuremice(
  imputationvariables,
  m = 25,
  maxit = 5,
  predictorMatrix = predM,
  n.core = 5,
  future.plan = "multisession",
  parallelseed = 123)

# =============================================================================
# GENERATE COMBINED TABLES
# =============================================================================

# Lithium
lithium_table_combined <- generate_combined_table(
  outcome = "lithium_prescribed",
  cc_data = model_data,
  imp_data = model_data_imp,
  vars_to_model = all_vars)

gtsave(lithium_table_combined,
       file = paste0(output_dir, "/lithium_regression_combined_", today(), ".docx"))

# Antidepressant
ad_table_combined <- generate_combined_table(
  outcome = "ad_only_prescribed",
  cc_data = model_data,
  imp_data = model_data_imp,
  vars_to_model = all_vars)

gtsave(ad_table_combined,
       file = paste0(output_dir, "/ad_regression_combined_", today(), ".docx"))

# =============================================================================
# FOREST PLOTS
# =============================================================================

# Create and save plots
# lithium_forest_plot <- generate_forest_plot(
#   outcome = "lithium_prescribed",
#   model_data_imp = model_data_imp,
#   vars_to_model = all_vars,
#   var_labels = var_labels,
#   plot_max = 1.3,
#   add_footnote = TRUE)
# 
# ad_forest_plot <- generate_forest_plot(
#   outcome = "ad_only_prescribed",
#   model_data_imp = model_data_imp,
#   vars_to_model = all_vars,
#   var_labels = var_labels,
#   plot_max = 1.5,
#   add_footnote = TRUE)
# 
# # Save as PNG and PDF
# for (plot_info in list(
#   list(plot = lithium_forest_plot, name = "lithium_forest_plot"),
#   list(plot = ad_forest_plot, name = "ad_forest_plot")
# )) {
#   ggsave(
#     filename = paste0(output_dir, "/", plot_info$name, "_", today(), ".png"),
#     plot = plot_info$plot,
#     width = 12, height = 10, dpi = 300
#   )
#   
#   pdf(
#     file = paste0(output_dir, "/", plot_info$name, "_", today(), ".pdf"),
#     width = 12, height = 10, bg = "white"
#   )
#   print(plot_info$plot)
#   dev.off()}

# =============================================================================
# Combined forest plot
# =============================================================================

# Reference headers and colors
REFERENCE_CHARS <- c("Sex", "IMD quintile", "Ethnicity", "Age at first diagnosis", "Diagnosis time-period")

CHARACTERISTIC_COLORS <- c(
  "Sex" = "#1f78b4",
  "IMD quintile" = "#33a02c",
  "Ethnicity" = "#e31a1c",
  "Age at first diagnosis" = "#ff7f00",
  "Diagnosis time-period" = "#6a3d9a")

# Assign colors function
assign_colors <- function(df) {
  df %>%
    mutate(char_group = case_when(
      characteristic %in% REFERENCE_CHARS ~ characteristic,
      TRUE ~ NA_character_
    )) %>%
    fill(char_group, .direction = "down") %>%
    mutate(
      color = CHARACTERISTIC_COLORS[char_group],
      color = ifelse(is_header, paste0(color, "B3"), color)
    ) %>%
    select(-char_group)}

# Full forest plot function #plot_title = "Treatment Prescribing Patterns"
plot_forest_faceted <- function(model_data_imp, vars_to_model, var_labels,
                                outcomes = c("lithium_prescribed", "ad_only_prescribed"),
                                digits_plot = 4) {
  
  # Header labels mapping
  HEADER_REF_LABELS <- c(
    "Sex" = "Sex (Reference: Female)",
    "IMD quintile" = "IMD quintile (Reference: 1 [Least deprived])",
    "Ethnicity" = "Ethnicity (Reference: White)",
    "Age at first diagnosis" = "Age at first diagnosis",
    "Diagnosis time-period" = "Diagnosis time-period (Reference: 2000–2005)"
  )
  
  all_outcome_results <- list()
  
  for (outcome in outcomes) {
    all_results <- list()
    
    for (var in vars_to_model) {
      label_list <- list()
      label_list[[var]] <- var_labels[[var]]
      
      mi_models <- fit_models(outcome, var, imp_data = model_data_imp)
      mi_tbl <- if (var == "time_period") mi_models$unadj else mi_models$adj
      
      tbl_df <- as.data.frame(
        tbl_regression(
          mi_tbl,
          exponentiate = TRUE,
          label = label_list,
          estimate_fun = purrr::partial(style_ratio, digits = digits_plot)
        )
      ) %>%
        clean_names() %>%
        separate(x95_percent_ci, into = c("conf_low", "conf_high"), sep = ", ") %>%
        mutate(
          across(c(or, conf_low, conf_high), as.numeric),
          is_header = characteristic %in% unname(var_labels),
          plot_label = ifelse(is_header, HEADER_REF_LABELS[characteristic], paste0("  ", characteristic)),
          ci_label = ifelse(!is.na(or),
                            paste0(format(round(or, 2), nsmall = 2),
                                   " (", format(round(conf_low, 2), nsmall = 2),
                                   "-", format(round(conf_high, 2), nsmall = 2), ")"),
                            "")
        ) %>%
        filter(!grepl("bipolar_year", characteristic)) %>%
        select(plot_label, ci_label, estimate = or, conf.low = conf_low,
               conf.high = conf_high, is_header, characteristic)
      
      all_results[[var]] <- tbl_df
    }
    
    mi_table <- do.call(rbind, all_results) %>%
      filter(!(is_header == FALSE & is.na(estimate))) %>%
      mutate(outcome = outcome)  # add outcome column
    
    all_outcome_results[[outcome]] <- mi_table
  }
  
  # Combine outcomes
  plot_data <- do.call(rbind, all_outcome_results) %>%
    assign_colors()
  
  plot_data <- plot_data %>%
    mutate(outcome = case_when(
      outcome == "lithium_prescribed" ~ "Lithium",
      outcome == "ad_only_prescribed" ~ "Antidepressant monotherapy",
      TRUE ~ as.character(outcome)
    )) %>%
    mutate(outcome = factor(outcome, levels = c("Lithium", "Antidepressant monotherapy")))
  
  # Step 2: Create SHARED row numbers (not per outcome)
  plot_data <- plot_data %>%
    mutate(char_group = ifelse(is_header, characteristic, NA_character_)) %>%
    fill(char_group, .direction = "down")
  
  plot_data <- plot_data %>%
    mutate(plot_label = gsub("\u2013", "-", plot_label))
  
  # Get row ordering from first outcome only
  first_outcome <- plot_data %>%
    filter(outcome == first(outcome)) %>%
    mutate(shared_row = rev(row_number())) %>%  # reverse for bottom-to-top
    select(char_group, characteristic, shared_row)
  
  # Apply same row numbers to all outcomes
  plot_data <- plot_data %>%
    left_join(first_outcome, by = c("char_group", "characteristic")) %>%
    mutate(show_label = outcome == "Lithium")
  
  check <- plot_data %>% 
    select(outcome, characteristic, shared_row) %>%
    pivot_wider(names_from = outcome, values_from = shared_row, names_prefix = "row_")
  print(check)
  
  # Step 3: Axis limits
  x_min <- min(plot_data$conf.low, na.rm = TRUE)
  x_max <- max(plot_data$conf.high, na.rm = TRUE)
  plot_min <- max(0, x_min - 0.08)
  plot_max <- max(x_max + 0.05, 1.18)
  arrow_threshold <- 1.05
  
  # Step 4: Force facet order: Lithium left, AD right
  plot_data$outcome <- factor(plot_data$outcome, levels = c("Lithium", "Antidepressant monotherapy"))
  
  # Step 5: Build the plot
  # Step 5: Build the plot
  p <- ggplot(plot_data, aes(y = shared_row)) +
    # Confidence intervals - normal (within both bounds)
    geom_segment(
      data = filter(plot_data, !is.na(estimate) & conf.high <= arrow_threshold & conf.low >= 0.31),
      aes(x = conf.low, xend = conf.high, yend = shared_row, color = color), 
      linewidth = 1
    ) +
    # Upper bound arrow only (right side)
    geom_segment(
      data = filter(plot_data, !is.na(estimate) & conf.high > arrow_threshold & conf.low >= 0.31),
      aes(x = conf.low, xend = arrow_threshold, yend = shared_row, color = color),
      linewidth = 1,
      arrow = arrow(angle = 20, length = unit(0.15, "cm"), type = "closed")
    ) +
    # Lower bound arrow only (left side)
    geom_segment(
      data = filter(plot_data, !is.na(estimate) & conf.low < 0.31 & conf.high <= arrow_threshold),
      aes(x = 0.31, xend = conf.high, yend = shared_row, color = color),
      linewidth = 1,
      arrow = arrow(angle = 20, length = unit(0.15, "cm"), type = "closed", ends = "first")
    ) +
    # Both bounds clipped (arrows on both ends)
    geom_segment(
      data = filter(plot_data, !is.na(estimate) & conf.low < 0.31 & conf.high > arrow_threshold),
      aes(x = 0.31, xend = arrow_threshold, yend = shared_row, color = color),
      linewidth = 1,
      arrow = arrow(angle = 20, length = unit(0.15, "cm"), type = "closed", ends = "both")
    ) +
    # Point estimates
    geom_point(
      data = filter(plot_data, !is.na(estimate)),
      aes(x = estimate, fill = color),
      color = "white", shape = 22, size = 4, stroke = 1.5
    ) +
    # Reference line
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.6) +
    # Left-hand labels (Lithium only)
    geom_text(
      data = filter(plot_data, show_label),
      aes(x = plot_min, y = shared_row, label = plot_label,
          fontface = ifelse(is_header, "bold","plain")),
      hjust = 0, size = 3.5, color = "black"
    ) +
    # Right CI labels
    geom_text(
      data = filter(plot_data, !is.na(estimate)),
      aes(x = plot_max, y = shared_row, label = ci_label),
      hjust = 1, size = 3.5
    ) +
    # Facet by outcome
    facet_wrap(~ outcome, ncol = 2, scales = "fixed") +
    scale_x_continuous(name = "Odds Ratio (95% CI)", limits = c(plot_min, plot_max), expand = c(0,0)) +
    scale_color_identity() +
    scale_fill_identity() +
    coord_cartesian(xlim = c(plot_min, plot_max), clip = "off") +
    theme_minimal(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "grey90", color = "grey60"),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray92", linewidth = 0.5),
      axis.line.x = element_line(color = "black", linewidth = 0.6),
      axis.ticks.x = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length.x = unit(0.15, "cm"),
      plot.margin = unit(c(0.25,1,1,1.5), "cm"),
      axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 10)),
      axis.text.x = element_text(size = 9, color = "black"),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 15)),
      plot.caption = element_text(size = 9, hjust = 0, color = "grey40", margin = margin(t = 10)),
      panel.spacing = unit(0.75, "lines"))
  #     ) +
  # #    ggtitle(plot_title) +
  #     labs(caption = "*All estimates are adjusted for diagnosis year and year², except time period (unadjusted due to collinearity).\nMissing values handled via multiple imputation (m=25).")
  
  print("Final ggplot object ready.")
  return(p)}

forest_faceted <- plot_forest_faceted(
  model_data_imp = model_data_imp,
  vars_to_model = all_vars,
  var_labels = var_labels,
  outcomes = c("lithium_prescribed", "ad_only_prescribed"))
#  plot_title = "Treatment Prescribing Patterns")

# 4. Save
ggsave(
  filename = paste0(output_dir, "/lithiumAD_forest_plot_", today(), ".png"),
  plot = forest_faceted,
  width = 45,    # cm
  height = 25,   # cm
  units = "cm",
  dpi = 300,
  bg = "white")

# 4. Save as PDF
ggsave(
  filename = paste0(output_dir, "/lithiumAD_forest_plot_", today(), ".pdf"),
  plot = forest_faceted,
  width = 45,    # cm
  height = 25,   # cm
  units = "cm",
  bg = "white")

## Overall risk and odds of each outcome
overall_risk_odds <- map_dfr(c("lithium_prescribed", "ad_only_prescribed"), ~ tibble(
  Outcome = .x,
  Sample_size = nrow(model_data),
  Risk = round(mean(model_data[[.x]], na.rm = TRUE) * 100, 1),
  Odds = round(exp(coef(glm(model_data[[.x]] ~ 1, family = "binomial"))[1]), 3)))
