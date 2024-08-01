# MODY-calculator-clinical-paper

This repository contains the code for

-   Make predictions for new patients.
-   Generating analysis and plots.

## Folder structure

1.   `data`: folder containing the data used for the analysis (empty as data cannot be shared)


2.   `figures`: folder containing figures used for the manuscript


3.   `model_development`: R code and model fitting objects (several models are included based on the different models fitted for the article)


4.   `model_predictions`: Vectors containing MODY probability estimations for all of the datasets used (`_full.rds` represents the full posterior dataset, its omission represents the summary results (2.5%, mean, 97.5%))


5.   `new_data_predictions`: R code + data + Excel document for making predictions in new datasets


## File explanation

1. `discrimination_curves.R`: this file contains the R code for calculating and plotting ROC and Precision-recall curves.