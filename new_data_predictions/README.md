# MODY probability predictions

> This contains the R code needed to make model predictions for new patients using R code.

The MODY calculator use a combination of two models to workout the probability of MODY.

An explanation on how to do predictions can be found below or in the file `prediction_functions.R`


Clinical information:


| **Field**            | **Requirements / Notes**                                                                     |
| -------------------- | -------------------------------------------------------------------------------------------- |
| **Current Age**      | - Must be between 1 and 120 years<br>- Must be equal to or greater than the Age of Diagnosis |
| **BMI (kg/m²)**      | - Use the most recent BMI<br>- Range: 10–70 kg/m²                                            |
| **HbA1c**            | - Use the most recent measurement<br>- 3–15% or 9–140 mmol/mol                               |
| **Age of Diagnosis** | - Must be between 1 and 35 years                                                             |




## Patients insulin-treated with 6-months of diagnosis

| **Test**                       | **Requirements / Notes**                                                                                              |
| ------------------------------ | --------------------------------------------------------------------------------------------------------------------- |
| **C-peptide Testing**          | - Use the result closest to calculator usage<br>- Plasma: 5–5000 pmol/L or 0.1–15 ng/ml<br>- UCPCR: 0.01–20 nmol/mmol |
| **Islet Autoantibody Testing** | - Use result closest to diagnosis<br>- If not tested near diagnosis, use first available result                       |



<details>
<summary>How to make predictions using R-software:</summary>
<br> 

To make predictions, you will need all of the files present on this folder [here](https://github.com/Exeter-Diabetes/MODY-calculator-clinical-paper/tree/main/new_data_predictions).

1.  Load the functions used for prediction.

``` r
# load functions
source("prediction_functions.R")
```

2.  Next, load the dataset containing the patient information.

``` r
# load dataset  
dataset <- ...
```

The dataset should be formatted in the following way:

| pardm<br>numeric | agerec<br>numeric | hba1c<br>numeric | agedx<br>numeric | sex<br>numeric | bmi<br>numeric | C<br>numeric | A<br>numeric |
|---------|---------|---------|---------|---------|---------|---------|---------|
| 1 - positive<br>0 - negative | \>1 or \<35 | \>3% or \< 15% | \>1 or \<120 | 1 - male<br>2 - female | \>14 or \<70 | 1 - positive<br>0 - negative | 1 - positive<br>0 - negative |

Note: 

- HbA1c - please include the most recent HbA1c measurement.
- BMI - please include the most recent BMI measurement.
- A - please include islet autoantibody results as close to diagnosis as possible. If not tested near to diagnosis, please use the first available islet autoantibody results as possible closest to diagnosis.


3.  Load the necessary model parameters.

``` r
# ## load posteriors
# rcs_parms <- readRDS("rcs_parms.rds")
# posterior_samples_T1D <- readRDS("type_1_model_posteriors_single_value.rds")
# 
# # ### create object to use for prediction
# posterior_samples_T1D_obj <- list(post = posterior_samples_T1D$samples)
# class(posterior_samples_T1D_obj) <- "T1D"
```

4.  Make predictions for new patients

``` r
## make predictions
posterior_predictions_T1D <- predict(posterior_samples_T1D_obj, dataset, rcs_parms) %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
  }) %>%
  bind_rows()
```

This code will produce a point prediction (`prob`), alongside a 95% credible interval (`LCI`-`UCI`).

<br>
</details>

## Patients non-/insulin-treated after 6-months of diagnosis

<details>
<summary>How to make predictions using R-software:</summary>
<br> 

To make predictions, you will need all of the files present on this folder [here](https://github.com/Exeter-Diabetes/MODY-calculator-clinical-paper/tree/main/new_data_predictions).

1.  Load the functions used for prediction.

``` r
# load functions
source("prediction_functions.R")
```

2.  Next, load the dataset containing the patient information.

``` r
# load dataset
dataset <- ...
```

The dataset should be formatted in the following way:

| pardm<br>numeric | agerec<br>numeric | hba1c<br>numeric | agedx<br>numeric | sex<br>numeric | bmi<br>numeric | insoroha<br>numeric |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| 1 - positive<br>0 - negative | \>1 or \<35 | \>3% or \<15% | \>1 or \<120 | 1 - male<br>2 - female | \>14 or \<70   | 1 - positive<br>0 - negative |

Note: 

- HbA1c - please include the most recent HbA1c measurement.
- BMI - please include the most recent BMI measurement.


3.  Load the necessary model parameters.

``` r
# ## load posteriors# 
# posterior_samples_T2D <- readRDS("type_2_model_posteriors_single_value.rds")
# 
# posterior_samples_T2D_obj <- list(post = posterior_samples_T2D$samples)
# class(posterior_samples_T2D_obj) <- "T2D"
```

4.  Make predictions for new patients

``` r
## make predictions
posterior_predictions_T2D <- predict(posterior_samples_T2D_obj, dataset) %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
  }) %>%
  bind_rows()
```

This code will produce a point prediction (`prob`), alongside a 95% credible interval (`LCI`-`UCI`).

<br>
</details>
