# MODY-calculator-clinical-paper

This repository contains the code for

-   Make predictions for new patients.
-   Generating analysis and plots.

## Folder structure

-   `data`: folder containing the data used for the analysis (empty as data cannot be shared)


-   `figures`: folder containing figures used for the manuscript


-   `model_development`: R code and model fitting objects (several models are included based on the different models fitted for the article)

  + `type_1_old_model.R` / `type_1_old_model_posteriors.rds`: R code and model posteriors for the old calculator model fitted to the case-control dataset in early-insulin-treated patients.
  + `type_2_old_model.R` / `type_2_old_model_posteriors.rds`: R code and model posteriors for the old calculator model fitted to the case-control dataset in non-early-insulin-treated patients.
  + `type_1_model.R` / `type_1_model_posteriors.rds`: R code and model posteriors for the new calculator model fitted to the UNITED dataset in early-insulin-treated patients (using the common MODY genes - HNF1a/HNF4a/GCK).
  + `type_2_model.R` / `type_2_model_posteriors.rds`: R code and model posteriors for the new calculator model fitted to the UNITED dataset in non-early-insulin-treated patients (using the common MODY genes - HNF1a/HNF4a/GCK).
  + `type_1_model_sensitivity_analysis.R` / `type_1_model_posteriors_sensitivity_analysis.rds`: R code and model posteriors for the new calculator model fitted to the UNITED dataset in early-insulin-treated patients (using the common MODY genes - HNF1a/HNF4a/GCK and only using antibody GAD information).
  + `rcs_parms.rds`: dataset of variable splines used for `type_1_model.R` and `type_1_model_sensitivity_analysis.R`.
  + `type_1_model_all_genes.R` / `type_1_model_posteriors_all_genes.rds`: R code and model posteriors for the new calculator model fitted to the UNITED dataset in early-insulin-treated patients (using the all MODY genes).
  + `rcs_parms_all_genes.rds`: dataset of variable splines used for `type_1_model_all_genes.R`.
  + `type_2_model_all_genes.R` / `type_2_model_posteriors_all_genes.rds`: R code and model posteriors for the new calculator model fitted to the UNITED dataset in non-early-insulin-treated patients (using the all MODY genes).


-   `model_predictions`: Vectors containing MODY probability estimations for all of the datasets used (`_full.rds` represents the full posterior dataset, its omission represents the summary results (2.5%, mean, 97.5%))


-   `new_data_predictions`: R code + data for making predictions in new datasets


