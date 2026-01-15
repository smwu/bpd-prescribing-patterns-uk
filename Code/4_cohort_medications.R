#Create a table listing  all antipsychotic, antidepressant, and mood stabiliser medication prescribed within the cohort between 2000 and 2022
# 
#
#
# Data:
#   - Data/bipolar_cohort_ads.rdata: Antidepressant prescriptions
#   - Data/bipolar_cohort_moodstabs.rdata: Mood stabiliser prescriptions
#   - Data/bipolar_cohort_aps.rdata: Antipsychotic prescriptions
#   - Data/bipolar_cohort.rdata: Patient information 
# Outputs:
#   - ingredients_formatted: formatted GT table listing all antipsychotic, antidepressant, and mood stabiliser medication prescribed within the cohort 
# ==============================================================================

# Clear memory
rm(list = ls())

library(lubridate)
library(gtsummary)
library(tidyverse)
library(ggforce)
library(gt)
library(officer)
library(flextable)

# Load data

load("Data/bipolar_cohort_ads.rdata")
load("Data/bipolar_cohort_moodstabs.rdata")
load("Data/bipolar_cohort_aps.rdata")
load("Data/bipolar_cohort.rdata")

#Filter our unknown gender
bipolar_cohort <- bipolar_cohort %>%
  filter(!gender %in% c("Indeterminate", "Unknown"))

# Combine valporate variants into "Valproate"
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
  mutate(drug_type = "ads")

aps <- bipolar_cohort_aps %>%
  select(patid, issuedate, ingredient) %>%
  mutate(drug_type = "aps")

moodstabs <- bipolar_cohort_moodstabs %>%
  select(patid, issuedate, ingredient) %>%
  mutate(drug_type = "moodstabs")

# Combine all prescriptions & remove duplicates
all_prescriptions <- bind_rows(ads, aps, moodstabs) %>%
  distinct()

# Replace codes with readable drug classes
ingredients_table <- all_prescriptions %>%
  mutate(drug_type = case_when(
    drug_type == "ads" ~ "Antidepressant",
    drug_type == "aps" ~ "Antipsychotic",
    drug_type == "moodstabs" ~ "Mood Stabiliser",
    TRUE ~ drug_type
  )) %>%
  select(ingredient, drug_type) %>%
  distinct()

# Wide format, keeping blanks instead of NA
ingredients_wide <- ingredients_table %>%
  arrange(drug_type, ingredient) %>%
  group_by(drug_type) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    names_from = drug_type,
    values_from = ingredient
  ) %>%
  select(Antidepressant, Antipsychotic, `Mood Stabiliser`) %>%
  mutate(across(everything(), ~replace_na(.x, "")))   # replace NA with blank


# Create gt table with grouped columns
ingredients_formatted <- ingredients_wide %>%
  gt() %>%
  tab_header(
    title = "List of Medications Prescribed to Patients with Bipolar Disorder in UK Primary Care: 2000 to 2022"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  )

#save as a .rtf file

gtsave(ingredients_formatted, "list_medications_gtable.rtf")

