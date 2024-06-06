#:--------------------------------------------------------
#   
# In this file we make predictions for each calculator for each dataset.
#
#:--------------------------------------------------------

# load libraries
library(nimble)
library(rms)
library(tidyverse)

# load functions needed for generating data
source("data/create_data.R")
source("new_data_predictions/prediction_functions.R")

### predict method for 'post' objects
predict.old_calculator_T1D <- function(object, newdata, ...) {
  
  ## check input objects
  stopifnot(class(object) == "old_calculator_T1D")
  stopifnot(is.data.frame(newdata))
  stopifnot(c("pardm", "agerec", "hba1c", "agedx", "sex") %in% colnames(newdata))
  x <- select(newdata, one_of(c("pardm", "agerec", "hba1c", "agedx", "sex")))
  stopifnot(all(map_chr(x, class) == "numeric"))
  
  ## convert posterior samples to matrix
  post <- as.matrix(object$post)
  
  ## set up data
  x <- as.matrix(x)
  x <- cbind(int = rep(1, nrow(x)), x)
  betas <- post[, match(c("beta0", paste0("beta[", 1:(ncol(x) - 1), "]")), colnames(post))]
  betas <- x %*% t(betas)
  
  ## do predictions ignoring tests
  preds <- exp(betas) / (1 + exp(betas))
  
  ## return posterior predictive samples
  preds <- t(preds)
  colnames(preds) <- NULL
  preds
}

### predict method for 'post' objects
predict.old_calculator_T2D <- function(object, newdata, ...) {
  
  ## check input objects
  stopifnot(class(object) == "old_calculator_T2D")
  stopifnot(is.data.frame(newdata))
  stopifnot(c("agedx", "bmi", "hba1c", "pardm", "agerec", "insoroha", "sex") %in% colnames(newdata))
  x <- select(newdata, one_of(c("agedx", "bmi", "hba1c", "pardm", "agerec", "insoroha", "sex")))
  stopifnot(all(map_chr(x, class) == "numeric"))
  
  ## convert posterior samples to matrix
  post <- as.matrix(object$post)
  
  ## set up data
  x <- as.matrix(x)
  x <- cbind(int = rep(1, nrow(x)), x)
  betas <- post[, match(c("beta0", paste0("beta[", 1:(ncol(x) - 1), "]")), colnames(post))]
  betas <- x %*% t(betas)
  
  ## do predictions ignoring tests
  preds <- exp(betas) / (1 + exp(betas))
  
  ## return posterior predictive samples
  preds <- t(preds)
  colnames(preds) <- NULL
  preds
}



# load datasets
## Load control-case dataset
dataset.case_control_type1 <- create_data(dataset = "case-control t1d")
dataset.case_control_type2 <- create_data(dataset = "case-control t2d")

## Load population representative dataset
dataset.UNITED_type1 <- create_data(dataset = "united t1d")
dataset.UNITED_type2 <- create_data(dataset = "united t2d")



## Load population representative dataset
dataset.UNITED_type1_gad <- create_data(dataset = "united t1d", biomarkers = "full") %>%
  
  # check if the antibody variable in question is recorded
  mutate(A = GAD) %>%
  
  # add biomarker latent variable
  mutate(T = ifelse(C == 0 | A == 1, 1, 0)) # T is 1 if Cn or Ap


dataset.UNITED_type1_gad_ia2 <- create_data(dataset = "united t1d", biomarkers = "full") %>%
  
  # check if the antibody variable in question is recorded
  mutate(
    A = ifelse(is.na(GAD) & is.na(IA2), NA,
               ifelse(!is.na(GAD) & (GAD == 1 & is.na(IA2)), 1,
                      ifelse(!is.na(GAD) & (GAD == 0 & is.na(IA2)), 0,
                             ifelse(!is.na(IA2) & (IA2 == 1 & is.na(GAD)), 1,
                                    ifelse(!is.na(IA2) & (IA2 == 0 & is.na(GAD)), 0,
                                           ifelse(IA2 == 1 | GAD == 1, 1, 0))))))
  ) %>%
  
  # add biomarker latent variable
  mutate(T = ifelse(C == 0 | A == 1, 1, 0)) # T is 1 if Cn or Ap




## Load referrals dataset

### load in Referral repository
current_path <- getwd()
setwd(dir = paste0(current_path, "/data/"))
#### download repository to get the latest updates
download.file(
  url = "https://github.com/Exeter-Diabetes/Julieanne-Pedro-MODY-Referrals/archive/master.zip",
  destfile = "meetingsR-master.zip"
)
#### unzip github repository
unzip(zipfile = "meetingsR-master.zip")
#### remove zipped file
file.remove("meetingsR-master.zip")
#### reverse working directory
setwd(current_path)
#### forget unnecessary variables
rm(current_path, create_data)

### load datasets
library(readxl)
source("data/Julieanne-Pedro-MODY-Referrals-main/data_formatting.R")

dataset.referral_type1 <- formatting(
  dataset = read_excel("data/mody_35yrs_08_02_2024.xlsx", guess_max = 1000000), 
  read_csv("data/mmoct11.csv"), 
  ethnicity_groups = read_excel("data/Julieanne-Pedro-MODY-Referrals-main/ethnicity_groups_ngs.xlsx", na = "null"), 
  ethnicity_groups_genetics = read_excel("data/Julieanne-Pedro-MODY-Referrals-main/ancestry_ngs_kp.xlsx", na = "null"), 
  ethnicity_labels_stated = read_excel("data/Julieanne-Pedro-MODY-Referrals-main/ethnicity_labels_stated.xlsx", na = "null"), 
  ethnicity_labels_genetics = read_excel("data/Julieanne-Pedro-MODY-Referrals-main/ethnicity_labels_genetics.xlsx", na = "null"), 
  biomarkers = "any", 
  diagnosis = FALSE, 
  type = "Type 1", 
  ethnicity = "White", 
  proband = "Proband", 
  gene_variable = FALSE
) %>%
  drop_na(c("sex", "bmi", "agedx", "hba1c", "pardm", "agerec", "bmi")) %>%   # what should be done about missing MODY testing?
  mutate(T = ifelse(C == 0 | A == 1, 1, 0)) # T is 1 if Cn or Ap

dataset.referral_type2 <- formatting(
  dataset = read_excel("data/mody_35yrs_08_02_2024.xlsx", guess_max = 1000000), 
  read_csv("data/mmoct11.csv"), 
  ethnicity_groups = read_excel("data/Julieanne-Pedro-MODY-Referrals-main/ethnicity_groups_ngs.xlsx", na = "null"), 
  ethnicity_groups_genetics = read_excel("data/Julieanne-Pedro-MODY-Referrals-main/ancestry_ngs_kp.xlsx", na = "null"), 
  ethnicity_labels_stated = read_excel("data/Julieanne-Pedro-MODY-Referrals-main/ethnicity_labels_stated.xlsx", na = "null"), 
  ethnicity_labels_genetics = read_excel("data/Julieanne-Pedro-MODY-Referrals-main/ethnicity_labels_genetics.xlsx", na = "null"), 
  diagnosis = FALSE, 
  type = "Type 2", 
  ethnicity = "White", 
  proband = "Proband", 
  gene_variable = FALSE
) %>%
  drop_na(c("sex", "bmi", "agedx", "insoroha", "hba1c", "pardm", "agerec"))    # what should be done about missing MODY testing?


#### remove downloaded folder
unlink("data/Julieanne-Pedro-MODY-Referrals-main", recursive = TRUE)

# load posteriors

## New calculator
rcs_parms <- readRDS("model_development/rcs_parms.rds")
posterior_samples_T1D <- readRDS("model_development/type_1_model_posteriors.rds")

# ### create object to use for prediction
posterior_samples_T1D_obj <- list(post = posterior_samples_T1D$samples)
class(posterior_samples_T1D_obj) <- "T1D"

posterior_samples_T2D <- readRDS("model_development/type_2_model_posteriors.rds")

posterior_samples_T2D_obj <- list(post = posterior_samples_T2D$samples)
class(posterior_samples_T2D_obj) <- "T2D"

### GAD info only model
posterior_samples_T1D_sensivity_analysis <- readRDS("model_development/type_1_model_posteriors_sensivity_analysis.rds")

# ### create object to use for prediction
posterior_samples_T1D_sensivity_analysis_obj <- list(post = posterior_samples_T1D_sensivity_analysis$samples)
class(posterior_samples_T1D_sensivity_analysis_obj) <- "T1D"

## Old calculator
posteriors_samples_old_T1D <- readRDS("model_development/type_1_old_model_posteriors.rds")

### create object to use for prediction
posteriors_samples_old_T1D <- list(post = posteriors_samples_old_T1D$samples)
class(posteriors_samples_old_T1D) <- "old_calculator_T1D"


posteriors_samples_old_T2D <- readRDS("model_development/type_2_old_model_posteriors.rds")

### create object to use for prediction
posteriors_samples_old_T2D <- list(post = posteriors_samples_old_T2D$samples)
class(posteriors_samples_old_T2D) <- "old_calculator_T2D"


## risk conversion (original method)
convert <- tibble(
  threshold = round(seq(0, 0.9, by = 0.1) * 100, 0), 
  PPVT1 = c(0.007, 0.019, 0.026, 0.040, 0.049, 0.064, 0.072, 0.082, 0.126, 0.494),
  PPVT2 = c(0.046, 0.151, 0.210, 0.244, 0.329, 0.358, 0.455, 0.580, 0.624, 0.755)
)



####################################################################################################################################################################################################################
####################################################################################################################################################################################################################

# Sensitivity analysis GAD only info

# Type 1 model

#UNITED all antibody
### Predictions from new model with T
interim <- as_tibble(as.matrix(select(dataset.UNITED_type1, pardm, agerec, hba1c, agedx, sex, bmi, C, A)))

predictions_dataset.UNITED_type1_sensitivity_analysis_with_T_full <- predict(posterior_samples_T1D_sensivity_analysis_obj, interim, rcs_parms)

#### save the predictions
saveRDS(predictions_dataset.UNITED_type1_sensitivity_analysis_with_T_full, "model_predictions/predictions_dataset.UNITED_type1_sensitivity_analysis_with_T_full.rds")

predictions_dataset.UNITED_type1_sensitivity_analysis_with_T <- predictions_dataset.UNITED_type1_sensitivity_analysis_with_T_full %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
  }) %>%
  bind_rows() 

#### save the predictions
saveRDS(predictions_dataset.UNITED_type1_sensitivity_analysis_with_T, "model_predictions/predictions_dataset.UNITED_type1_sensitivity_analysis_with_T.rds")


# UNITED gad dataset
### Predictions from new model with T
interim <- as_tibble(as.matrix(select(dataset.UNITED_type1_gad, pardm, agerec, hba1c, agedx, sex, bmi, C, A)))

predictions_dataset.UNITED_type1_gad_sensitivity_analysis_with_T_full <- predict(posterior_samples_T1D_sensivity_analysis_obj, interim, rcs_parms)

#### save the predictions
saveRDS(predictions_dataset.UNITED_type1_gad_sensitivity_analysis_with_T_full, "model_predictions/predictions_dataset.UNITED_type1_gad_sensitivity_analysis_with_T_full.rds")

predictions_dataset.UNITED_type1_gad_sensitivity_analysis_with_T <- predictions_dataset.UNITED_type1_gad_sensitivity_analysis_with_T_full %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
  }) %>%
  bind_rows() 

#### save the predictions
saveRDS(predictions_dataset.UNITED_type1_gad_sensitivity_analysis_with_T, "model_predictions/predictions_dataset.UNITED_type1_gad_sensitivity_analysis_with_T.rds")

####################################################################################################################################################################################################################
####################################################################################################################################################################################################################


# Type 1 model

### Predictions from new model no T
interim <- as_tibble(as.matrix(select(dataset.case_control_type1 %>%
                                        mutate(C = NA, A = NA), pardm, agerec, hba1c, agedx, sex, bmi, C, A)))


predictions_dataset.case_control_type1_no_T_full <- predict(posterior_samples_T1D_obj, interim, rcs_parms)

#### save the predictions
saveRDS(predictions_dataset.case_control_type1_no_T_full, "model_predictions/predictions_dataset.case_control_type1_no_T_full.rds")

predictions_dataset.case_control_type1_no_T <- predictions_dataset.case_control_type1_no_T_full %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
  }) %>%
  bind_rows() 

#### save the predictions
saveRDS(predictions_dataset.case_control_type1_no_T, "model_predictions/predictions_dataset.case_control_type1_no_T.rds")


#:-------------



interim <- as_tibble(as.matrix(select(dataset.referral_type1 %>%
                                        mutate(C = NA, A = NA), pardm, agerec, hba1c, agedx, sex, bmi, C, A)))

predictions_dataset.referral_type1_no_T_full <- predict(posterior_samples_T1D_obj, interim, rcs_parms)

#### save the predictions
saveRDS(predictions_dataset.referral_type1_no_T_full, "model_predictions/predictions_dataset.referral_type1_no_T_full.rds")

predictions_dataset.referral_type1_no_T <- predictions_dataset.referral_type1_no_T_full %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
  }) %>%
  bind_rows() 

#### save the predictions
saveRDS(predictions_dataset.referral_type1_no_T, "model_predictions/predictions_dataset.referral_type1_no_T.rds")



#:-------------



interim <- as_tibble(as.matrix(select(dataset.UNITED_type1 %>%
                                        mutate(C = NA, A = NA), pardm, agerec, hba1c, agedx, sex, bmi, C, A)))

predictions_dataset.UNITED_type1_no_T_full <- predict(posterior_samples_T1D_obj, interim, rcs_parms) 

#### save the predictions
saveRDS(predictions_dataset.UNITED_type1_no_T_full, "model_predictions/predictions_dataset.UNITED_type1_no_T_full.rds")

predictions_dataset.UNITED_type1_no_T <- predictions_dataset.UNITED_type1_no_T_full %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
  }) %>%
  bind_rows() 

#### save the predictions
saveRDS(predictions_dataset.UNITED_type1_no_T, "model_predictions/predictions_dataset.UNITED_type1_no_T.rds")



#:-------------



### Predictions from new model with T

interim <- as_tibble(as.matrix(select(dataset.referral_type1, pardm, agerec, hba1c, agedx, sex, bmi, C, A)))

predictions_dataset.referral_type1_with_T_full <- predict(posterior_samples_T1D_obj, interim, rcs_parms)

#### save the predictions
saveRDS(predictions_dataset.referral_type1_with_T_full, "model_predictions/predictions_dataset.referral_type1_with_T_full.rds")

predictions_dataset.referral_type1_with_T <- predictions_dataset.referral_type1_with_T_full %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
  }) %>%
  bind_rows() 

#### save the predictions
saveRDS(predictions_dataset.referral_type1_with_T, "model_predictions/predictions_dataset.referral_type1_with_T.rds")




#:-------------



interim <- as_tibble(as.matrix(select(dataset.UNITED_type1, pardm, agerec, hba1c, agedx, sex, bmi, C, A)))

predictions_dataset.UNITED_type1_with_T_full <- predict(posterior_samples_T1D_obj, interim, rcs_parms) 

#### save the predictions
saveRDS(predictions_dataset.UNITED_type1_with_T_full, "model_predictions/predictions_dataset.UNITED_type1_with_T_full.rds")

predictions_dataset.UNITED_type1_with_T <- predictions_dataset.UNITED_type1_with_T_full %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
  }) %>%
  bind_rows() 

#### save the predictions
saveRDS(predictions_dataset.UNITED_type1_with_T, "model_predictions/predictions_dataset.UNITED_type1_with_T.rds")



## GAD

interim <- as_tibble(as.matrix(select(dataset.UNITED_type1_gad, pardm, agerec, hba1c, agedx, sex, bmi, C, A)))

predictions_dataset.UNITED_type1_gad_with_T_full <- predict(posterior_samples_T1D_obj, interim, rcs_parms) 

#### save the predictions
saveRDS(predictions_dataset.UNITED_type1_gad_with_T_full, "model_predictions/predictions_dataset.UNITED_type1_gad_with_T_full.rds")

predictions_dataset.UNITED_type1_gad_with_T <- predictions_dataset.UNITED_type1_gad_with_T_full %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
  }) %>%
  bind_rows() 

#### save the predictions
saveRDS(predictions_dataset.UNITED_type1_gad_with_T, "model_predictions/predictions_dataset.UNITED_type1_gad_with_T.rds")


#:-------------


## GAD & IA2

interim <- as_tibble(as.matrix(select(dataset.UNITED_type1_gad_ia2, pardm, agerec, hba1c, agedx, sex, bmi, C, A)))

predictions_dataset.UNITED_type1_gad_ia2_with_T_full <- predict(posterior_samples_T1D_obj, interim, rcs_parms) 

#### save the predictions
saveRDS(predictions_dataset.UNITED_type1_gad_ia2_with_T_full, "model_predictions/predictions_dataset.UNITED_type1_gad_ia2_with_T_full.rds")

predictions_dataset.UNITED_type1_gad_ia2_with_T <- predictions_dataset.UNITED_type1_gad_ia2_with_T_full %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
  }) %>%
  bind_rows() 

#### save the predictions
saveRDS(predictions_dataset.UNITED_type1_gad_ia2_with_T, "model_predictions/predictions_dataset.UNITED_type1_gad_ia2_with_T.rds")



#:-------------


### Predictions from old model

interim <- as_tibble(as.matrix(select(dataset.case_control_type1, pardm, agerec, hba1c, agedx, sex)))

predictions_dataset.case_control_type1_old_full <- predict(posteriors_samples_old_T1D, interim) 

#### save the predictions
saveRDS(predictions_dataset.case_control_type1_old_full, "model_predictions/predictions_dataset.case_control_type1_old_full.rds")

predictions_dataset.case_control_type1_old <- predictions_dataset.case_control_type1_old_full%>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x))
  }) %>%
  bind_rows()

for (i in 1:nrow(predictions_dataset.case_control_type1_old)) {
  prob <- predictions_dataset.case_control_type1_old[i,]
  prob <- 10 * ((round(prob  * 100, 0)) %/% 10)
  if (prob == 100) {
    predictions_dataset.case_control_type1_old[i,] <- convert$PPVT1[10]
  } else {
    predictions_dataset.case_control_type1_old[i,] <- convert$PPVT1[convert$threshold == prob]
  }
}

#### save the predictions
saveRDS(predictions_dataset.case_control_type1_old, "model_predictions/predictions_dataset.case_control_type1_old.rds")



#:-------------




interim <- as_tibble(as.matrix(select(dataset.referral_type1, pardm, agerec, hba1c, agedx, sex)))

predictions_dataset.referral_type1_old_full <- predict(posteriors_samples_old_T1D, interim) 

#### save the predictions
saveRDS(predictions_dataset.referral_type1_old_full, "model_predictions/predictions_dataset.referral_type1_old_full.rds")

predictions_dataset.referral_type1_old <- predictions_dataset.referral_type1_old_full%>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x))
  }) %>%
  bind_rows()

for (i in 1:nrow(predictions_dataset.referral_type1_old)) {
  prob <- predictions_dataset.referral_type1_old[i,]
  prob <- 10 * ((round(prob  * 100, 0)) %/% 10)
  if (prob == 100) {
    predictions_dataset.referral_type1_old[i,] <- convert$PPVT1[10]
  } else {
    predictions_dataset.referral_type1_old[i,] <- convert$PPVT1[convert$threshold == prob]
  }
}

#### save the predictions
saveRDS(predictions_dataset.referral_type1_old, "model_predictions/predictions_dataset.referral_type1_old.rds")




#:-------------



interim <- as_tibble(as.matrix(select(dataset.UNITED_type1, pardm, agerec, hba1c, agedx, sex)))

predictions_dataset.UNITED_type1_old_full <- predict(posteriors_samples_old_T1D, interim) 

#### save the predictions
saveRDS(predictions_dataset.UNITED_type1_old_full, "model_predictions/predictions_dataset.UNITED_type1_old_full.rds")

predictions_dataset.UNITED_type1_old <- predictions_dataset.UNITED_type1_old_full %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x))
  }) %>%
  bind_rows()

for (i in 1:nrow(predictions_dataset.UNITED_type1_old)) {
  prob <- predictions_dataset.UNITED_type1_old[i,]
  prob <- 10 * ((round(prob  * 100, 0)) %/% 10)
  if (prob == 100) {
    predictions_dataset.UNITED_type1_old[i,] <- convert$PPVT1[10]
  } else {
    predictions_dataset.UNITED_type1_old[i,] <- convert$PPVT1[convert$threshold == prob]
  }
}

#### save the predictions
saveRDS(predictions_dataset.UNITED_type1_old, "model_predictions/predictions_dataset.UNITED_type1_old.rds")



#:-------------




# Type 2 model

### Predictions from new model

interim <- as_tibble(as.matrix(select(dataset.case_control_type2, pardm, agerec, hba1c, agedx, sex, bmi, insoroha)))

predictions_dataset.case_control_type2_new_full <- predict(posterior_samples_T2D_obj, interim) 

#### save the predictions
saveRDS(predictions_dataset.case_control_type2_new_full, "model_predictions/predictions_dataset.case_control_type2_new_full.rds")

predictions_dataset.case_control_type2_new <- predictions_dataset.case_control_type2_new_full%>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
  }) %>%
  bind_rows() 

#### save the predictions
saveRDS(predictions_dataset.case_control_type2_new, "model_predictions/predictions_dataset.case_control_type2_new.rds")



#:-------------




interim <- as_tibble(as.matrix(select(dataset.referral_type2, pardm, agerec, hba1c, agedx, sex, bmi, insoroha)))

predictions_dataset.referral_type2_new_full <- predict(posterior_samples_T2D_obj, interim) 

#### save the predictions
saveRDS(predictions_dataset.referral_type2_new_full, "model_predictions/predictions_dataset.referral_type2_new_full.rds")

predictions_dataset.referral_type2_new <- predictions_dataset.referral_type2_new_full %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
  }) %>%
  bind_rows() 

#### save the predictions
saveRDS(predictions_dataset.referral_type2_new, "model_predictions/predictions_dataset.referral_type2_new.rds")



#:-------------



interim <- as_tibble(as.matrix(select(dataset.UNITED_type2, pardm, agerec, hba1c, agedx, sex, bmi, insoroha)))

predictions_dataset.UNITED_type2_new_full <- predict(posterior_samples_T2D_obj, interim) 

#### save the predictions
saveRDS(predictions_dataset.UNITED_type2_new_full, "model_predictions/predictions_dataset.UNITED_type2_new_full.rds")

predictions_dataset.UNITED_type2_new <- predictions_dataset.UNITED_type2_new_full %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
  }) %>%
  bind_rows() 

#### save the predictions
saveRDS(predictions_dataset.UNITED_type2_new, "model_predictions/predictions_dataset.UNITED_type2_new.rds")



#:-------------



### Predictions from old model

interim <- as_tibble(as.matrix(select(dataset.case_control_type2, agedx, bmi, hba1c, pardm, agerec, insoroha, sex)))

predictions_dataset.case_control_type2_old_full <- predict(posteriors_samples_old_T2D, interim) 

#### save the predictions
saveRDS(predictions_dataset.case_control_type2_old_full, "model_predictions/predictions_dataset.case_control_type2_old_full.rds")

predictions_dataset.case_control_type2_old <- predictions_dataset.case_control_type2_old_full %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x))
  }) %>%
  bind_rows()

for (i in 1:nrow(predictions_dataset.case_control_type2_old)) {
  prob <- predictions_dataset.case_control_type2_old[i,]
  prob <- 10 * ((round(prob  * 100, 0)) %/% 10)
  if (prob == 100) {
    predictions_dataset.case_control_type2_old[i,] <- convert$PPVT2[10]
  } else {
    predictions_dataset.case_control_type2_old[i,] <- convert$PPVT2[convert$threshold == prob]
  }
}

#### save the predictions
saveRDS(predictions_dataset.case_control_type2_old, "model_predictions/predictions_dataset.case_control_type2_old.rds")



#:-------------





interim <- as_tibble(as.matrix(select(dataset.referral_type2, agedx, bmi, hba1c, pardm, agerec, insoroha, sex)))

predictions_dataset.referral_type2_old_full <- predict(posteriors_samples_old_T2D, interim) 

#### save the predictions
saveRDS(predictions_dataset.referral_type2_old_full, "model_predictions/predictions_dataset.referral_type2_old_full.rds")

predictions_dataset.referral_type2_old <- predictions_dataset.referral_type2_old_full%>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x))
  }) %>%
  bind_rows()

for (i in 1:nrow(predictions_dataset.referral_type2_old)) {
  prob <- predictions_dataset.referral_type2_old[i,]
  prob <- 10 * ((round(prob  * 100, 0)) %/% 10)
  if (prob == 100) {
    predictions_dataset.referral_type2_old[i,] <- convert$PPVT2[10]
  } else {
    predictions_dataset.referral_type2_old[i,] <- convert$PPVT2[convert$threshold == prob]
  }
}

#### save the predictions
saveRDS(predictions_dataset.referral_type2_old, "model_predictions/predictions_dataset.referral_type2_old.rds")




#:-------------



interim <- as_tibble(as.matrix(select(dataset.UNITED_type2, agedx, bmi, hba1c, pardm, agerec, insoroha, sex)))

predictions_dataset.UNITED_type2_old_full <- predict(posteriors_samples_old_T2D, interim) 

#### save the predictions
saveRDS(predictions_dataset.UNITED_type2_old_full, "model_predictions/predictions_dataset.UNITED_type2_old_full.rds")

predictions_dataset.UNITED_type2_old <- predictions_dataset.UNITED_type2_old_full%>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x))
  }) %>%
  bind_rows()

for (i in 1:nrow(predictions_dataset.UNITED_type2_old)) {
  prob <- predictions_dataset.UNITED_type2_old[i,]
  prob <- 10 * ((round(prob  * 100, 0)) %/% 10)
  if (prob == 100) {
    predictions_dataset.UNITED_type2_old[i,] <- convert$PPVT2[10]
  } else {
    predictions_dataset.UNITED_type2_old[i,] <- convert$PPVT2[convert$threshold == prob]
  }
}

#### save the predictions
saveRDS(predictions_dataset.UNITED_type2_old, "model_predictions/predictions_dataset.UNITED_type2_old.rds")






