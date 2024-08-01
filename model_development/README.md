# Model development folder:

This folder contains R code and model fitting objects (several models are included based on the different models fitted for the article).


1. `type_1_old_model.R` / `type_1_old_model_posteriors.rds`: R code and model posteriors for the old calculator model fitted to the case-control dataset in early-insulin-treated patients.
2. `type_2_old_model.R` / `type_2_old_model_posteriors.rds`: R code and model posteriors for the old calculator model fitted to the case-control dataset in non-early-insulin-treated patients.
3. `type_1_model.R` / `type_1_model_posteriors.rds`: R code and model posteriors for the new calculator model fitted to the UNITED dataset in early-insulin-treated patients (using the common MODY genes - HNF1a/HNF4a/GCK).
4. `type_2_model.R` / `type_2_model_posteriors.rds`: R code and model posteriors for the new calculator model fitted to the UNITED dataset in non-early-insulin-treated patients (using the common MODY genes - HNF1a/HNF4a/GCK).
5. `type_1_model_sensitivity_analysis.R` / `type_1_model_posteriors_sensitivity_analysis.rds`: R code and model posteriors for the new calculator model fitted to the UNITED dataset in early-insulin-treated patients (using the common MODY genes - HNF1a/HNF4a/GCK and only using antibody GAD information).
6. `rcs_parms.rds`: dataset of variable splines used for `type_1_model.R` and `type_1_model_sensitivity_analysis.R`.
7. `type_1_model_all_genes.R` / `type_1_model_posteriors_all_genes.rds`: R code and model posteriors for the new calculator model fitted to the UNITED dataset in early-insulin-treated patients (using the all MODY genes).
8. `rcs_parms_all_genes.rds`: dataset of variable splines used for `type_1_model_all_genes.R`.
9. `type_2_model_all_genes.R` / `type_2_model_posteriors_all_genes.rds`: R code and model posteriors for the new calculator model fitted to the UNITED dataset in non-early-insulin-treated patients (using the all MODY genes).
