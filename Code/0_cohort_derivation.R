# -----------------------
# BIPOLAR DISORDER PRESCRIBING TRENDS 2000-2022
# -----------------------
# CPRD 2023 data
# Last run: 29/10/2025
# 
# DERIVE COHORT
# 
# Data: 
#   - Data/PatSMI_C.Rdata: Bipolar diagnosis date data
#   - Data/SMIcohort.Rdata: SMI patient cohort information
#   - Data/antipsychotics_combined.Rdata: Antipsychotic prescription preprocessed data
#   - Data/SMImoodstabilisers_C.Rdata: Mood stabiliser prescription preprocessed data
#   - Data/SMIantidepressants.Rdata: Antidepressant prescription preprocessed data
# 
# Outputs:
#   - Data/bipolar_cohort.rdata: Bipolar disorder cohort with patient information
#   - Data/bipolar_cohort_smidiagnoses.rdata: Diagnosis date information for bipolar cohort
#   - Data/bipolar_cohort_aps.rdata: Antipsychotic prescriptions for bipolar cohort
#   - Data/bipolar_cohort_moodstabs.rdata: Mood stabiliser prescriptions for bipolar cohort
#   - Data/bipolar_cohort_ads.rdata: Antidepressant prescriptions for bipolar cohort
# ------------------------

# Clear memory
rm(list = ls())

# Packages
invisible(sapply(c("dplyr", "tidyr", "stringr", "forcats", "lubridate", "readr", "data.table", 
                   "purrr", "rlang", "reshape2", "ggplot2", "gtsummary", "gt", "tidylog"), 
                 library, character.only = TRUE))

# DATA MANAGEMENT ####

# Identify first bipolar disorder diagnosis date
load("Data/PatSMI_C.Rdata")

bipolar_first <- PatSMI_C %>%
  filter(group == "bipolar") %>%
  mutate(readterm = if_else(readterm == "SMI" & medcode == "53840",
                            "[X]Other bipolar affective disorders", readterm)) %>%
  filter(!str_detect(readterm, regex("lithium|remission", ignore_case = TRUE))) %>%
  group_by(patid) %>%
  summarise(first_bipolar_date = min(eventdate), .groups = "drop") %>%
  left_join(PatSMI_C, by = c("patid", "first_bipolar_date" = "eventdate")) %>%
  mutate(episode_type = case_when( # See if we can extract info on episode type
    str_detect(readterm, regex("hypomania|hypomanic", ignore_case = TRUE)) ~ "hypomania",
    str_detect(readterm, regex("\\bmania\\b|\\bmanic\\b(?![- ]depressive)", ignore_case = TRUE)) ~ "mania",
    str_detect(readterm, regex("mixed", ignore_case = TRUE)) ~ "mixed/rapid cycling",
    str_detect(readterm, regex("type II|Bipolar 2|bipolar II", ignore_case = TRUE)) ~ "type 2",
    str_detect(readterm, regex("type I|Bipolar I|bipolar I", ignore_case = TRUE)) ~ "type 1",
    str_detect(readterm, regex("(?<!manic[- ])depressed|(?<!manic[- ])depression|(?<!manic[- ])depressive|(?<!manic[- ])depressn|(?<!manic[- ])severe depres|(?<!manic[- ])sev depress", ignore_case = TRUE)) ~ "depressed",
    TRUE ~ "not specified")) %>%
  group_by(patid) %>%
  summarise(
    first_bipolar_date = first(first_bipolar_date),
    types_all = paste(sort(unique(episode_type[episode_type != "not specified"])), collapse = ", "),
    types_all = if_else(types_all == "", "not specified", types_all),
    types_single = case_when(
      types_all == "not specified" ~ "not specified",
      str_detect(types_all, "mixed") ~ "mixed",
      str_detect(types_all, "type 1") ~ "type 1",
      str_detect(types_all, "type 2") ~ "type 2",
      str_detect(types_all, "mania") & str_detect(types_all, "depress") ~ "mixed",
      str_detect(types_all, "hypomania") & str_detect(types_all, "depress") ~ "mixed",
      str_detect(types_all, "\\bmania\\b") & str_detect(types_all, "hypomania") ~ "mania",
      TRUE ~ types_all),
    .groups = "drop")

# Cohort
load("Data/SMIcohort.Rdata")

# =========================
# 1. Data processing
# =========================

bipolar_cohort <- SMIcohort %>%
  select(patid, pracid, gender, yob, regstartdate, regenddate, lcd, source, region, deathdate, ethnicity_cat_cprdhes, 
         patprac_2019imd_quintile, pat_2019imd_quintile, prac_2019imd_quintile) %>%
  left_join(bipolar_first, by = "patid") %>%
  mutate(bipolar_year = year(first_bipolar_date),
             age_at_first_bipolar = bipolar_year - yob) %>%
  mutate(pre_followup_days  = as.numeric(first_bipolar_date - regstartdate),
         post_followup_days = as.numeric(pmin(deathdate, regenddate, lcd, na.rm = TRUE) - first_bipolar_date),
         pre_followup_years  = pre_followup_days / 365.25,
         post_followup_years = post_followup_days / 365.25,
         pre_cat_1yr  = ifelse(pre_followup_days < 365, "<1 yr", "≥1 yr"),
         pre_cat_2yr  = ifelse(pre_followup_days < 730, "<2 yrs", "≥2 yrs"),
         post_cat_1yr = ifelse(!is.na(deathdate) & (deathdate - first_bipolar_date) <= 365, "Died ≤1 yr", "≥1 yr"),
         post_cat_2yr = ifelse(!is.na(deathdate) & (deathdate - first_bipolar_date) <= 730, "Died ≤2 yrs", "≥2 yrs")) %>%
  # Remove missing bipolar diagnosis
  filter(!is.na(first_bipolar_date)) %>%
  # Not censored before first bipolar date
  filter(post_followup_days >= 0) %>%
  # BP diagnosis between 2000/01/01 and <2023/01/01
  filter(first_bipolar_date >= '2000-01-01', first_bipolar_date < '2023-01-01') %>%
  # Bipolar diagnosis between 18-99 years of age
  filter(age_at_first_bipolar >= 18, age_at_first_bipolar < 100) %>%
  # At least 365 days followup before diagnosis
  filter(pre_followup_days >= 365) %>%
  # At least 365 days followup after diagnosis OR died within 365 days after diagnosis
  filter(post_followup_days >= 365 | (!is.na(deathdate) & (deathdate - first_bipolar_date) <= 365))

# SMI diagnoses
bipolar_cohort_smidiagnoses <- bipolar_cohort %>%
  select(patid) %>%
  inner_join(PatSMI_C, by = "patid") %>%
  filter(!grepl("(?i)lithium", readterm)) %>%
  filter(eventdate < '2023-01-01') %>%
  distinct()

rm(PatSMI_C)

# Antipsychotics
load("Data/antipsychotics_combined.Rdata")

bipolar_cohort_aps <- bipolar_cohort %>%
  select(patid) %>%
  inner_join(antipsychotics_combined, by = "patid") %>%
  filter(issuedate < '2023-01-01') %>%
  filter(AP != "Prochlorperazine") %>%
  select(-estnhscost, -probobsid)

rm(antipsychotics_combined)

table(bipolar_cohort_aps$AP)
n_distinct(bipolar_cohort_aps$patid)

# Mood stabilisers
load("Data/SMImoodstabilisers_C.Rdata")

bipolar_cohort_moodstabs <- bipolar_cohort %>%
  select(patid) %>%
  inner_join(SMImoodstabilisers_C, by = "patid") %>%
  filter(issuedate < '2023-01-01') %>%
  select(-estnhscost, -probobsid)

# Combine valproate variants into "Valproate"
bipolar_cohort_moodstabs <- bipolar_cohort_moodstabs %>%
  mutate(ingredient = case_when(ingredient %in% c("Valproate", "Valproate/Valproic", "Valproic") ~ "Valproate",
                                TRUE ~ ingredient),
         MoodStabiliser = case_when(MoodStabiliser %in% c("Sodium valproate", "Valproic acid", "Valproate/Valproic") ~ "Valproate",
                                TRUE ~ MoodStabiliser))
         
rm(SMImoodstabilisers_C)

table(bipolar_cohort_moodstabs$MoodStabiliser)
n_distinct(bipolar_cohort_moodstabs$patid)

# Antidepressants
load("Data/SMIantidepressants.Rdata")

bipolar_cohort_ads <- bipolar_cohort %>%
  select(patid) %>%
  inner_join(SMIantidepressants, by = "patid") %>%
  filter(issuedate < '2023-01-01') %>%
  select(-estnhscost, -probobsid)

rm(SMIantidepressants)

table(bipolar_cohort_ads$Antidepressant)
n_distinct(bipolar_cohort_ads$patid)

gc()

# Recode gender
bipolar_cohort <- bipolar_cohort %>%
  mutate(gender = case_when(gender == "Male"   ~ "Male",
                            gender == "Female" ~ "Female",
                            TRUE ~ NA_character_))

# Save data files
save(bipolar_cohort, file = "Data/bipolar_cohort.rdata")
save(bipolar_cohort_smidiagnoses, file = "Data/bipolar_cohort_smidiagnoses.rdata")
save(bipolar_cohort_aps, file = "Data/bipolar_cohort_aps.rdata")
save(bipolar_cohort_moodstabs, file = "Data/bipolar_cohort_moodstabs.rdata")
save(bipolar_cohort_ads, file = "Data/bipolar_cohort_ads.rdata")