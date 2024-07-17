# MODY probability predictions

> This contains the R code needed to make model predictions for new patients using R code.

The MODY calculator use a combination of two models to workout the probability of MODY.

An explanation on how to do predictions can be found below or in the file `prediction_functions.R`

## Patients insulin-treated with 6-months of diagnosis

<details>
<summary>How to make predictions using R-software:</summary>
<br> 

To make predictions, you will need all of the files present on this folder [here](https://github.com/Exeter-Diabetes/MODY-calculator-clinical-paper/tree/main/new_data_predictions).

1.  Load the functions used for prediction.

``` r
# load functions
source("prediction_functions.R")
```

2.  Next, load the data containing the patient information.

``` r
# load dataset  
data <- ...
```

The data should be formatted in the following way:

| pardm<br>numeric | agerec<br>numeric | hba1c<br>numeric | agedx<br>numeric | sex<br>numeric | bmi<br>numeric | C<br>numeric | A<br>numeric |
|---------|---------|---------|---------|---------|---------|---------|---------|
| 1 or 0 | \>1 or \<35 | \>3% or \< 15% | \>1 or \<120 | 1 - male<br>2 - female | \>14 or \<70 | 1 or 0 | 1 or 0 |

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

2.  Next, load the data containing the patient information.

``` r
# load dataset
data <- ...
```

The data should be formatted in the following way:

| pardm<br>numeric | agerec<br>numeric | hba1c<br>numeric | agedx<br>numeric | sex<br>numeric | bmi<br>numeric | insoroha<br>numeric |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| 1 or 0 | \>1 or \<35 | \>3% or \<15% | \>1 or \<120 | 1 - male<br>2 - female | \>14 or \<70   | 1 or 0 |

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
