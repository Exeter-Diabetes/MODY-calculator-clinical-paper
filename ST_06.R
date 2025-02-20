##########################################################################################
#Supplementary Table 6

#This code produces Supplementary Table 6

#These tables contain the diagnostic & discriminatory criteria resulting 
#from using either a 5, 10, 15, 20, 25, or 30% testing threshold
#applied to the external validation cohort used in this paper for both early- and not-
#early-insulin-treated participants

###############################################################################################

#SET WORKING DIRECTORY -----------------------------------------------------------------------------------
setwd("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper")

#load libraries
library(tidyverse)
library(writexl)

#load functions
source("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper/new_data_predictions/prediction_functions.R")
source("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper/data/create_data.R")

#load data
# MyDiabetes ------------------------------------------------------------------------------------
#load T1D
load("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper/MY_T1D.RData")
#load T2D 
load("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper/MY_T2D.RData")
# LIMIT TO ONLY WHITE ETHNICITY
MY_T1D <- MY_T1D %>%
  filter(ethnicity == "White")
MY_T1D <- MY_T1D %>%
  rename(
    hba1c_mmol = hba1c
  )

MY_T1D <- MY_T1D %>%
  rename(
    hba1c = hba1c_perc, 
    bmi = BMI
  )

#T2D data
MY_T2D <- MY_T2D %>%
  filter(ethnicity == "White")

MY_T2D <- MY_T2D %>%
  rename(
    hba1c_mmol = hba1c
  )

MY_T2D <- MY_T2D %>%
  rename(
    hba1c = hba1c_perc, 
    bmi = BMI
  )

# UNITED paediatric ---------------------------------------------------------------------
dataset.UNITED_p <- create_data(dataset = "united t1d pediatrics", commonmody = FALSE)

#need to change M=NA to M=0
dataset.UNITED_p <- dataset.UNITED_p %>%
  mutate(M = ifelse(is.na(M) == TRUE, 0, M)) %>%
  filter(!is.na(T))
#checked if worked: should have M=1 (n=7) & M=0 (n=1164)
table(dataset.UNITED_p$M)

#T2D
dataset.UNITED_type2p <- dataset.UNITED_p %>%
  filter(tti == "" | tti == "Greater than 12 months")
table(dataset.UNITED_type2p$M)

#T1D
dataset.UNITED_type1p <- dataset.UNITED_p %>%
  filter(!(id %in% dataset.UNITED_type2p$id))
table(dataset.UNITED_type1p$M)

# External datasets -----------------------------------------------------------------------

# Early-insulin-treated
dataset.MYDIABETES_type1 <- MY_T1D %>%
  mutate(study = "MYDIABETES") %>%
  select(MY_ID, study, agerec, agedx, sex, bmi, pardm, insoroha, hba1c, C, A, M, Gene, biomark_status) 

dataset.UNITED_type1p <- dataset.UNITED_type1p %>%
  mutate(study = "UNITED paediatric") %>%
  filter(!is.na(T))
dataset_type1 <- full_join(dataset.MYDIABETES_type1, dataset.UNITED_type1p, by = c("study","agerec", "agedx", "sex", "bmi", "pardm", "insoroha", "hba1c", "C", "A", "M", "Gene")) %>%
  mutate(M = ifelse(is.na(M), 0, M), 
         Gene = ifelse(Gene == "", NA, Gene), 
         biomark_status = ifelse(is.na(biomark_status), 
                                 ifelse(C == 1 & A == 0, 0, 1),
                                 biomark_status)) 
# Not-early-insulin-treated
dataset.MYDIABETES_type2 <- MY_T2D %>%
  mutate(study = "MYDIABETES") %>%
  select(MY_ID, study, agerec, agedx, sex, bmi, pardm, insoroha, hba1c, C, A, M, Gene, biomark_status) 

dataset.UNITED_type2p <- dataset.UNITED_type2p %>%
  mutate(study = "UNITED paediatric") %>%
  select(id, study, agerec, agedx, sex, bmi, pardm, insoroha, hba1c, C, A, M, Gene) 

dataset_type2 <- full_join(dataset.MYDIABETES_type2, dataset.UNITED_type2p, by = c("study","agerec", "agedx", "sex", "bmi", "pardm", "insoroha", "hba1c", "C", "A", "M", "Gene")) %>%
  mutate(M = ifelse(is.na(M), 0, M), 
         Gene = ifelse(Gene == "", NA, Gene), 
         biomark_status = ifelse(is.na(biomark_status), 
                                 ifelse(C == 1 & A == 0, 0, 1),
                                 biomark_status)) 

#2. APPLY MODELS
setwd("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper/new_data_predictions")

# Early-insulin-treated ------------------------------------------------------------------------------------------
# ## load posteriors
rcs_parms <- readRDS("rcs_parms.rds")

##T1D/EARLY-INSULIN-TREATED MODEL 
posterior_samples_T1D <- readRDS("type_1_model_posteriors_thin_100.rds")
#
# ### create object to use for prediction
posterior_samples_T1D_obj <- list(post = posterior_samples_T1D$samples)
class(posterior_samples_T1D_obj) <- "T1D"

## make predictions
#this yields simplified dataset with 3 columns: 
#prob (which is the mean probability across the bayesian probability iterations), 
#lower confidence interval and upper confidence interval
posterior_predictions_T1D <- predict(posterior_samples_T1D_obj, dataset_type1, rcs_parms) %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, prob = 0.025), UCI = quantile(x, prob = 0.975))
  }) %>%
  bind_rows()

#posterior_predictions_external_val_t1d_full <- predict(posterior_samples_T1D_obj, dataset_type1, rcs_parms)
#save(posterior_predictions_external_val_t1d_full, file = "posterior_predictions_external_val_t1d_full.rds")

#MERGE TO dataset_type1
dataset_type1 <- cbind(dataset_type1, posterior_predictions_T1D)

# Not-early-insulin-treated  -----------------------------------------------------------------------------------------
posterior_samples_T2D <- readRDS("type_2_model_posteriors_thin_100.rds")
# 
posterior_samples_T2D_obj <- list(post = posterior_samples_T2D$samples)
class(posterior_samples_T2D_obj) <- "T2D"

## make predictions
posterior_predictions_T2D <- predict(posterior_samples_T2D_obj, dataset_type2) %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, prob = 0.025), UCI = quantile(x, prob = 0.975))
  }) %>%
  bind_rows()

#this yields the full (huge) probability estimations and iterations - takes approx an hour to run (with nothing else functioning - not the case)
#only ever do once
#posterior_predictions_external_val_t2d_full <- predict(posterior_samples_T2D_obj, dataset_type2)
#save(posterior_predictions_external_val_t2d_full, file = "posterior_predictions_external_val_t2d_full.rds")


#MERGE TO dataset_type2
dataset_type2 <- cbind(dataset_type2, posterior_predictions_T2D)

#Joint tables 
External_joint <- full_join(dataset_type1, dataset_type2)

# Save as table --------------------------------------------------------------------------------------------

setwd("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper")
# Early-insulin-treated
### 5%
dataset_5PERC <- External_joint %>%
  summarise(totalover = sum(prob >= 0.05),
            ncasespickedup = sum(prob >= 0.05 & M ==1),
            PPV = (sum(prob >= 0.05 & M ==1)/sum(prob >= 0.05))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.05 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.05 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.05 & M == 0)/sum(prob < 0.05))*100,
            Sensitivity = (sum(prob >= 0.05 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.05 & M == 0)/sum(M==0))*100)

### 10%
dataset_10PERC <- External_joint %>%
  summarise(totalover = sum(prob >= 0.1),
            ncasespickedup = sum(prob >= 0.1 & M ==1),
            PPV = (sum(prob >= 0.1 & M ==1)/sum(prob >= 0.1))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.1 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.1 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.1 & M == 0)/sum(prob < 0.1))*100,
            Sensitivity = (sum(prob >= 0.1 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.1 & M == 0)/sum(M==0))*100)


### 15%
dataset_15PERC <- External_joint %>%
  summarise(totalover = sum(prob >= 0.15),
            ncasespickedup = sum(prob >= 0.15 & M ==1),
            PPV = (sum(prob >= 0.15 & M ==1)/sum(prob >= 0.15))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.15 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.15 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.15 & M == 0)/sum(prob < 0.15))*100,
            Sensitivity = (sum(prob >= 0.15 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.15 & M == 0)/sum(M==0))*100)

### 20%
dataset_20PERC <- External_joint %>%
  summarise(totalover = sum(prob >= 0.2),
            ncasespickedup = sum(prob >= 0.2 & M ==1),
            PPV = (sum(prob >= 0.2 & M ==1)/sum(prob >= 0.2))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.2 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.2 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.2 & M == 0)/sum(prob < 0.2))*100,
            Sensitivity = (sum(prob >= 0.2 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.2 & M == 0)/sum(M==0))*100)

### 25%
dataset_25PERC <- External_joint %>%
  summarise(totalover = sum(prob >= 0.25),
            ncasespickedup = sum(prob >= 0.25 & M ==1),
            PPV = (sum(prob >= 0.25 & M ==1)/sum(prob >= 0.25))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.25 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.25 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.25 & M == 0)/sum(prob < 0.25))*100,
            Sensitivity = (sum(prob >= 0.25 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.25 & M == 0)/sum(M==0))*100)

### 30%
dataset_30PERC <- External_joint %>%
  summarise(totalover = sum(prob >= 0.3),
            ncasespickedup = sum(prob >= 0.3 & M ==1),
            PPV = (sum(prob >= 0.3 & M ==1)/sum(prob >= 0.3))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.3 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.3 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.3 & M == 0)/sum(prob < 0.3))*100,
            Sensitivity = (sum(prob >= 0.3 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.3 & M == 0)/sum(M==0))*100)

dataset_THRESHOLDS <- rbind(dataset_5PERC, 
                            dataset_10PERC,
                            dataset_15PERC, 
                            dataset_20PERC, 
                            dataset_25PERC, 
                            dataset_30PERC)
write_xlsx(dataset_THRESHOLDS,"Supp_Table6.xlsx")


