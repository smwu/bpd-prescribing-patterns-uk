# -----------------------
# BIPOLAR DISORDER PRESCRIBING TRENDS 2000-2022
# -----------------------
# CPRD 2023 data
# 
# Baseline characteristics ####
# Last updated: 20/10/2025
# 
# Data:
#   - Data/bipolar_cohort.rdata: Patient information 
#   - Data/SMIModSevDepression_C.Rdata: Pre-existing depression information for patients
#   - Data/SMIDiabetes_C.Rdata: Pre-existing diabetes data
#   - Data/SMIHypertension_C.Rdata: Pre-existing hypertension data
#   - Data/SMILiverDisease_C.Rdata: Pre-existing liver disease data
#   - Data/SMIRenalDisease_C.Rdata: Pre-existing renal disease data
#   - Data/SMICharlson_C.Rdata": Charlson comorbidity index
#   - Data/SMIElixhauser_C.Rdata: Elixhauser comorbidity index, used to assess alcohol and substance misuse at baseline
# 
# Outputs:
#   - Outputs/baseline_characteristics: Table of baseline characteristics 
#   - Outputs/followup_summary: Summary table of follow-up time
#   - Outputs/diagnosis_year_histogram: Histogram of diagnosis year
# -----------------------


# Clear memory
rm(list = ls())

# Libraries
library(dplyr)
library(gtsummary)
library(gt)
library(tidyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(scales)


# Set output directory
output_dir <- "Outputs/"
wd <- getwd()

# Load cohort
load("Data/bipolar_cohort.rdata")

# DATA MANAGEMENT ####

# Pre-existing depression
load("Data/SMIModSevDepression_C.Rdata")

depression_flags <- bipolar_cohort %>%
  inner_join(SMIModSevDepression_C, by = "patid") %>%
  group_by(patid, first_bipolar_date) %>%
  summarise(
    depression_ever_before_bipolar = as.integer(any(eventdate < first_bipolar_date)),
    depression_inlast1y_before_bipolar = as.integer(any(
      eventdate < first_bipolar_date & eventdate >= first_bipolar_date - 365)),
    .groups = "drop")

bipolar_cohort <- bipolar_cohort %>%
  left_join(depression_flags, by = c("patid", "first_bipolar_date")) %>%
  mutate(
    depression_ever_before_bipolar = ifelse(is.na(depression_ever_before_bipolar), 0, depression_ever_before_bipolar),
    depression_inlast1y_before_bipolar = ifelse(is.na(depression_inlast1y_before_bipolar), 0, depression_inlast1y_before_bipolar))

rm(depression_flags, SMIModSevDepression_C)

# Other comorbidities
extract_baseline_comorbidities <- function(diagnosis_data, var_name) {
  
  # Construct the file name to load the diagnosis data
  load(paste0("Data/", diagnosis_data, ".Rdata"))
  
  # CPRD data
  diagnosis_prev <- bipolar_cohort %>%
    select(patid, first_bipolar_date) %>%
    inner_join(get(diagnosis_data), by = "patid") %>%
    filter(difftime(eventdate, first_bipolar_date, units = "days") < 0) %>%
    distinct(patid) %>%
    mutate(!!var_name := 1)
  
  remove(diagnosis_data)
  gc()
  
  return(diagnosis_prev)}

# Call the function for baseline diagnoses
diabetes_prev <- extract_baseline_comorbidities("SMIDiabetes_C", "priordiabetes")
htn_prev <- extract_baseline_comorbidities("SMIHypertension_C", "priorhypertension")
liver_prev <- extract_baseline_comorbidities("SMILiverDisease_C", "priorliverdisease")
renal_prev <- extract_baseline_comorbidities("SMIRenalDisease_C", "priorrenaldisease")

load("Data/SMICharlson_C.Rdata")
load("Data/SMIElixhauser_C.Rdata")

# Substance abuse at baseline
substance_prev <- bipolar_cohort %>%
  select(patid, first_bipolar_date) %>%
  inner_join(SMIElixhauser_C %>%
               filter(condition == "Drug abuse" & !grepl("(?i)alcohol misuse", readterm)), by = "patid") %>%
  filter(difftime(eventdate, first_bipolar_date, units = "days") < 0) %>%
  select(patid) %>%
  distinct() %>%
  mutate(priorsubstanceabuse = 1)

#Alcohol abuse at baseline
alcohol_prev <- bipolar_cohort %>%
  select(patid, first_bipolar_date) %>%
  inner_join(SMIElixhauser_C %>%
               filter(condition == "Alcohol abuse"),
             by = "patid") %>%
  filter(difftime(eventdate, first_bipolar_date, units = "days") < 0) %>%
  select(patid) %>%
  distinct() %>%
  mutate(prioralcoholabuse = 1)

remove(SMICharlson_C, SMIElixhauser_C)
gc()

df_list <- list(bipolar_cohort, 
                substance_prev, alcohol_prev, htn_prev, diabetes_prev, liver_prev, renal_prev)

bipolar_cohort <- df_list %>% purrr::reduce(full_join, by = "patid") %>%
  mutate(across(30:35, ~ if_else(is.na(.), "No", if_else(. == 1, "Yes", "No")))) # set any missing concomitant medication values to 'No' and convert to factor

rm(df_list, substance_prev, alcohol_prev, htn_prev, diabetes_prev, liver_prev, renal_prev)

# TABLES AND FIGURES ####

# Baseline characteristics table
missing_pat_imd <- sum(bipolar_cohort$pat_2019imd_quintile %in% NA & !is.na(bipolar_cohort$prac_2019imd_quintile))

# Define comorbidity variable labels for consistent formatting
comorbidity_labels <- c("Moderate/severe depression \\(ever pre-existing\\)", 
                        "Moderate/severe depression \\(in last 1y\\)",
                        "Alcohol misuse", "Substance misuse", 
                        "Diabetes", "Hypertension", 
                        "Liver disease", "Renal disease")

chars_table <- bipolar_cohort %>%
  mutate(ethnicity_cat_cprdhes = str_to_title(ethnicity_cat_cprdhes),
         types_single = str_to_title(types_single)) %>%
  select(age_at_first_bipolar, gender, ethnicity_cat_cprdhes, region, patprac_2019imd_quintile, 
         bipolar_year, types_single, depression_ever_before_bipolar, depression_inlast1y_before_bipolar, 
         prioralcoholabuse, priorsubstanceabuse, priordiabetes, priorhypertension, priorliverdisease, priorrenaldisease) %>%
  tbl_summary(
    label = list(age_at_first_bipolar = "Age at first bipolar disorder diagnosis, median (IQR)",
                 gender = "Sex, n (%)",
                 ethnicity_cat_cprdhes = "Ethnicity, n (%)",
                 bipolar_year = "Year of first bipolar disorder diagnosis, median (IQR)",
                 types_single = "Episode type, n (%)",
                 region = "Region, n (%)",
                 patprac_2019imd_quintile = "IMD Quintile, n (%)",
                 depression_ever_before_bipolar = "Moderate/severe depression (ever pre-existing)",
                 depression_inlast1y_before_bipolar = "Moderate/severe depression (in last 1y)",
                 priordiabetes = "Diabetes",
                 priorhypertension = "Hypertension",
                 priorliverdisease = "Liver disease", 
                 priorrenaldisease = "Renal disease",
                 prioralcoholabuse = "Alcohol misuse", 
                 priorsubstanceabuse = "Substance misuse"),
    digits = list(all_categorical() ~ c(0, 1))) %>%
  # Remove comma formatting from year values
  modify_table_body(
    ~.x %>%
      mutate(stat_0 = if_else(
        variable == "bipolar_year",
        str_remove_all(stat_0, ","),
        stat_0
      ))
  ) %>%
  # Convert to tibble to add Comorbidities header row
  as_tibble() %>%
  {
    depression_row <- which(.[[1]] == "Moderate/severe depression (ever pre-existing)")
    
    # Insert "Comorbidities" heading before first comorbidity
    if(length(depression_row) > 0) {
      bind_rows(
        slice(., 1:(depression_row - 1)),
        tibble(!!names(.)[1] := "Comorbidities, n (%)", !!names(.)[2] := ""),
        slice(., depression_row:n())
      )
    } else {
      .
    }
  } %>%
  # Remove asterisks from column names
  rename_with(~str_remove_all(., "\\*\\*"), everything()) %>%
  # Convert to gt table
  gt() %>%
  tab_style(style = cell_text(weight = "bold"), # format spanners bold
            locations = cells_column_labels(everything())) %>%
# Replace NA values with empty strings
  sub_missing(missing_text = "") %>%
  # Indent all category rows (those that don't contain "n (%)" or "median (IQR)" and aren't "Comorbidities")
  text_transform(
    locations = cells_body(
      columns = Characteristic,
      rows = !grepl("n \\(%\\)|median \\(IQR\\)|Comorbidities, n (%)", Characteristic)
    ),
    fn = function(x) {
      paste0("\u00A0\u00A0\u00A0\u00A0", x)  # 4 non-breaking spaces
    }
  ) %>%
  # Bold column headers
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(columns = everything())
  ) %>%
  # Bold variable labels (excluding individual comorbidities)
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = Characteristic,
      rows = (grepl("n \\(%\\)|median \\(IQR\\)", Characteristic) & 
                !grepl(paste(comorbidity_labels, collapse = "|"), Characteristic)) |
        Characteristic == "Comorbidities, n (%)"
    )
  ) %>%
  # Footnote for IMD variable
  tab_footnote(
    footnote = paste0("Quintile of the 2019 English Index of Multiple Deprivation - defined according to the patient's postcode or, where this was missing (n=", missing_pat_imd, "), the primary care practice postcode."),
    locations = cells_body(columns = Characteristic, rows = Characteristic == "IMD Quintile, n (%)")
  ) %>%
  # Footnote for Comorbidities section
  tab_footnote(
    footnote = "Comorbidities determined at any point in the patient's medical history prior to the index date.",
    locations = cells_body(columns = Characteristic, rows = Characteristic == "Comorbidities, n (%)"))
  
# Save
gtsave(chars_table,
       file = paste0(output_dir, "/baseline_characteristics_", today(), ".docx"))

rm(missing_pat_imd, comorbidity_labels, chars_table)

# Follow-up time summary table

# Helper function to format summary statistics
format_summary_stats <- function(x) {
  q25 <- formatC(quantile(x, 0.25, na.rm = TRUE), format = "f", digits = 1)
  q75 <- formatC(quantile(x, 0.75, na.rm = TRUE), format = "f", digits = 1)
  
  paste0(
    formatC(mean(x, na.rm = TRUE), format = "f", digits = 1), " ± ",
    formatC(sd(x, na.rm = TRUE), format = "f", digits = 1), "; ",
    formatC(median(x, na.rm = TRUE), format = "f", digits = 1), " (",
    q25, "–", q75, ")")}

# Build the summary table
summary_table <- bind_rows(
  
  # Section 1: Cohort size and deaths
  bipolar_cohort %>%
    summarise(
      Section = "Cohort size & deaths",
      Measure = "Total N",
      Value = format(n(), big.mark = ",")
    ),
  
  # Section 2: Pre-diagnosis follow-up (validation of ≥1 yr criterion)
  bipolar_cohort %>%
    summarise(
      Section = "Pre-diagnosis follow-up",
      Measure = "Pre-follow-up (yrs)",
      Value = format_summary_stats(pre_followup_years)
    ),
  
  # Pre-diagnosis 1-year category
  bipolar_cohort %>%
    count(pre_cat_1yr) %>%
    mutate(
      Section = "Pre-diagnosis follow-up",
      Measure = paste0("  ", pre_cat_1yr),
      Value = paste0(format(n, big.mark = ","), " (", round(100 * n / sum(n), 1), "%)")
    ) %>%
    select(Section, Measure, Value),
  
  # Pre-diagnosis ≥2 yrs category only, percentage relative to full cohort
  bipolar_cohort %>%
    count(pre_cat_2yr) %>%
    filter(pre_cat_2yr == "≥2 yrs") %>%
    mutate(
      Section = "Pre-diagnosis follow-up",
      Measure = paste0("  ", pre_cat_2yr),
      Value = paste0(
        format(n, big.mark = ","), " (", 
        round(100 * n / nrow(bipolar_cohort), 1), "%)")) %>%
    select(Section, Measure, Value),
  
  # Section 3: Post-diagnosis follow-up
  bipolar_cohort %>%
    summarise(
      Section = "Post-diagnosis follow-up",
      Measure = "Post-follow-up (yrs)",
      Value = format_summary_stats(post_followup_years)),
  
  # Post-diagnosis 1-year categories
  bipolar_cohort %>%
    count(post_cat_1yr) %>%
    mutate(
      Section = "Post-diagnosis follow-up",
      Measure = paste0("  ", post_cat_1yr),
      Value = paste0(format(n, big.mark = ","), " (", round(100 * n / sum(n), 1), "%)")
    ) %>%
    select(Section, Measure, Value),
  
  # Post-diagnosis 2-year categories
  bipolar_cohort %>%
    count(post_cat_2yr) %>%
    mutate(
      Section = "Post-diagnosis follow-up",
      Measure = paste0("  ", post_cat_2yr),
      Value = paste0(format(n, big.mark = ","), " (", round(100 * n / sum(n), 1), "%)")
    ) %>%
    select(Section, Measure, Value))

# Render the gt table
summary_table <- summary_table %>%
  gt(groupname_col = "Section") %>%
  tab_header(title = "Follow-up Summary") %>%
  cols_label(
    Measure = "",
    Value = "") %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()) %>%
  tab_options(
    row_group.background.color = "#F5F5F5",
    row_group.font.weight = "bold")

# Save
gtsave(summary_table,
       file = paste0(output_dir, "/followup_summary_", today(), ".docx"))

rm(summary_table)

# Year of diagnosis plot
diagnosis_year <- ggplot(bipolar_cohort, aes(x = bipolar_year)) +
  geom_histogram(
    breaks = seq(2000, 2022, by = 1),
    color = "black",
    fill = "#1f77b4"
  ) +
  scale_x_continuous(
    breaks = seq(2000, 2022, by = 2),
    limits = c(2000, 2023),  # ensure all years are included
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    labels = comma,  # add thousand separators
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Year of first recorded diagnosis",
    y = "Number of patients"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),  # rotate x-axis labels
    panel.grid.minor = element_blank())

ggsave(filename = paste0("diagnosis_year_histogram_", today(), ".pdf"),
       path = output_dir, plot = diagnosis_year,
       width = 8, height = 5, units = "in")

rm(diagnosis_year)
