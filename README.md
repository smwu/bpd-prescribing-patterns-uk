# bpd-prescribing-patterns-uk

## Psychiatric prescribing patterns in patients with newly diagnosed bipolar disorder in the United Kingdom: 2000 to 2022

This repository contains code for analyzing patterns of antidepressant, antipsychotic, and mood stabiliser prescribing among people newly diagnosed with bipolar disorder in UK primary care. The data source is [Clinical Practice Research Datalink (CPRD)](https://www.cprd.com/).

### Code

The following R scripts were used to perform data processing and analyses:

-   `0_cohort_derivation.R`: Derive and process data for cohort of patients with first bipolar diagnosis between 2000 and 2022.
-   `1_baseline_characteristics.R`: Create table of baseline characteristics for patient cohort.
-   `2_monthly_line_plot.R`: Create line plots of the proportion of patients prescribed each medication class in the 12 months before and after first-recorded BPD diagnosis.
-   `3_sankey_plot.R`: Create Sankey plots and heatmaps of Patterns of treatment persistence and medication class switching within one year of first-recorded BPD diagnosis.
-   `4_cohort_medications.R`: Create a summary table listing all antipsychotics, antidepressants, and mood stabilisers prescibed to the cohort between 2000 and 2022.
-   `5_regression_models.R`: Run regression models examining associations with lithium and antidepressant monotherapy prescription post-diagnosis, and save output.

### Code lists

The `Code_Lists` folder contains CPRD GOLD and Aurum code lists that were used to define antidepressant, antipsychotic, and mood stabiliser prescriptions.

### Shiny app

The `bpd-prescribing-uk-shiny` folder contains code (`app.R`) and aggregated, de-identified data used to create an interactive Shiny web application that was developed to visualise treatment trajectories over the first year following diagnosis.

### Outputs

The `Outputs` folder contains saved analysis outputs that are created from the code.
