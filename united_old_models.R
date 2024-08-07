## UNITED_OLD_MODELS -------------------------------------------------------------------

# load libraries ---------------------------------------------------------------------------------
library(tidyverse)
library(writexl)

## load functions needed for generating data ------------------------------------
source("data/create_data.R")

#load UNITED data ---------------------------------------------------------------------------------
## 3 genes ----------------------------------------------------------------------------------------------------
## Load model predictions
predictions_dataset.UNITED_type1_old <- readRDS("model_predictions/predictions_dataset.UNITED_type1_old.rds")
predictions_dataset.UNITED_type2_old <- readRDS("model_predictions/predictions_dataset.UNITED_type2_old.rds")

##load cohort data
dataset.UNITED_type1_old <- create_data(dataset = "united t1d") %>%
  
  ## if MODY testing missing, change to 0
  mutate(M = ifelse(is.na(M), 0, M))

dataset.UNITED_type2_old <- create_data(dataset = "united t2d")

##merge cohort and predictions data
UNITED_type1_old <- cbind(dataset.UNITED_type1_old,predictions_dataset.UNITED_type1_old)
UNITED_type2_old <- cbind(dataset.UNITED_type2_old,predictions_dataset.UNITED_type2_old)

## all genes -----------------------------------------------------------------------------------------------
## Load model predictions
predictions_dataset.UNITED_type1_all_genes_old <- readRDS("model_predictions/predictions_dataset.UNITED_type1_all_genes_old.rds")
predictions_dataset.UNITED_type2_all_genes_old <- readRDS("model_predictions/predictions_dataset.UNITED_type2_all_genes_old.rds")

## load cohort data
dataset.UNITED_type1_all_genes <- create_data(dataset = "united t1d", commonmody = FALSE) %>%
  ## if MODY testing missing, change to 0
  mutate(M = ifelse(is.na(M), 0, M))

dataset.UNITED_type2_all_genes <- create_data(dataset = "united t2d", commonmody = FALSE)

##merge cohort and predictions data
UNITED_type1_all_genes_old <- cbind(dataset.UNITED_type1_all_genes,predictions_dataset.UNITED_type1_all_genes_old)
UNITED_type2_all_genes_old <- cbind(dataset.UNITED_type2_all_genes,predictions_dataset.UNITED_type2_all_genes_old)




##Thresholds -----------------------------------------------------------------------------
### 3 genes --------------------------------------------------------------------------------------
#### T1D ------------------------------------------------------------------------------------
##### 5%
UNITED_type1_old_5PERC <- UNITED_type1_old %>%
  summarise(totalover = sum(prob >= 0.05),
            ncasespickedup = sum(prob >= 0.05 & M ==1),
            PPV = (sum(prob >= 0.05 & M ==1)/sum(prob >= 0.05))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.05 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.05 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.05 & M == 0)/sum(prob < 0.05))*100,
            Sensitivity = (sum(prob >= 0.05 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.05 & M == 0)/sum(M==0))*100)
##### 10%
UNITED_type1_old_10PERC <- UNITED_type1_old %>%
  summarise(totalover = sum(prob >= 0.1),
            ncasespickedup = sum(prob >= 0.1 & M ==1),
            PPV = (sum(prob >= 0.1 & M ==1)/sum(prob >= 0.1))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.1 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.1 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.1 & M == 0)/sum(prob < 0.1))*100,
            Sensitivity = (sum(prob >= 0.1 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.1 & M == 0)/sum(M==0))*100)

### 20%
UNITED_type1_old_20PERC <- UNITED_type1_old %>%
  summarise(totalover = sum(prob >= 0.2),
            ncasespickedup = sum(prob >= 0.2 & M ==1),
            PPV = (sum(prob >= 0.2 & M ==1)/sum(prob >= 0.2))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.2 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.2 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.2 & M == 0)/sum(prob < 0.2))*100,
            Sensitivity = (sum(prob >= 0.2 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.2 & M == 0)/sum(M==0))*100)

UNITED_type1_old_THRESHOLDS <- rbind(UNITED_type1_old_5PERC, UNITED_type1_old_10PERC, UNITED_type1_old_20PERC)
write_xlsx(UNITED_type1_old_THRESHOLDS,"UNITED_type1_old_THRESHOLDS_table.xlsx")

##Thresholds -----------------------------------------------------------------------------
##T2D ----------------------------------------------------------------------------------------
### 5%
UNITED_type2_old_5PERC <- UNITED_type2_old %>%
  summarise(totalover = sum(prob >= 0.05),
            ncasespickedup = sum(prob >= 0.05 & M ==1),
            PPV = (sum(prob >= 0.05 & M ==1)/sum(prob >= 0.05))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.05 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.05 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.05 & M == 0)/sum(prob < 0.05))*100,
            Sensitivity = (sum(prob >= 0.05 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.05 & M == 0)/sum(M==0))*100)
### 10%
UNITED_type2_old_10PERC <- UNITED_type2_old %>%
  summarise(totalover = sum(prob >= 0.1),
            ncasespickedup = sum(prob >= 0.1 & M ==1),
            PPV = (sum(prob >= 0.1 & M ==1)/sum(prob >= 0.1))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.1 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.1 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.1 & M == 0)/sum(prob < 0.1))*100,
            Sensitivity = (sum(prob >= 0.1 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.1 & M == 0)/sum(M==0))*100)

### 20%
UNITED_type2_old_20PERC <- UNITED_type2_old %>%
  summarise(totalover = sum(prob >= 0.2),
            ncasespickedup = sum(prob >= 0.2 & M ==1),
            PPV = (sum(prob >= 0.2 & M ==1)/sum(prob >= 0.2))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.2 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.2 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.2 & M == 0)/sum(prob < 0.2))*100,
            Sensitivity = (sum(prob >= 0.2 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.2 & M == 0)/sum(M==0))*100)

UNITED_type2_old_THRESHOLDS <- rbind(UNITED_type2_old_5PERC, UNITED_type2_old_10PERC, UNITED_type2_old_20PERC)
write_xlsx(UNITED_type2_old_THRESHOLDS,"UNITED_type2_old_THRESHOLDS_table.xlsx")

### all genes ----------------------------------------------------------------------------------------
#### T1D ------------------------------------------------------------------------------------
### 5%
UNITED_type1_all_genes_old_5PERC <- UNITED_type1_all_genes_old %>%
  summarise(totalover = sum(prob >= 0.05),
            ncasespickedup = sum(prob >= 0.05 & M ==1),
            PPV = (sum(prob >= 0.05 & M ==1)/sum(prob >= 0.05))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.05 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.05 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.05 & M == 0)/sum(prob < 0.05))*100,
            Sensitivity = (sum(prob >= 0.05 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.05 & M == 0)/sum(M==0))*100)
### 10%
UNITED_type1_all_genes_old_10PERC <- UNITED_type1_all_genes_old %>%
  summarise(totalover = sum(prob >= 0.1),
            ncasespickedup = sum(prob >= 0.1 & M ==1),
            PPV = (sum(prob >= 0.1 & M ==1)/sum(prob >= 0.1))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.1 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.1 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.1 & M == 0)/sum(prob < 0.1))*100,
            Sensitivity = (sum(prob >= 0.1 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.1 & M == 0)/sum(M==0))*100)

### 20%
UNITED_type1_all_genes_old_20PERC <- UNITED_type1_all_genes_old %>%
  summarise(totalover = sum(prob >= 0.2),
            ncasespickedup = sum(prob >= 0.2 & M ==1),
            PPV = (sum(prob >= 0.2 & M ==1)/sum(prob >= 0.2))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.2 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.2 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.2 & M == 0)/sum(prob < 0.2))*100,
            Sensitivity = (sum(prob >= 0.2 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.2 & M == 0)/sum(M==0))*100)

UNITED_type1_all_genes_old_THRESHOLDS <- rbind(UNITED_type1_all_genes_old_5PERC, UNITED_type1_all_genes_old_10PERC, UNITED_type1_all_genes_old_20PERC)
write_xlsx(UNITED_type1_all_genes_old_THRESHOLDS,"UNITED_type1_all_genes_old_THRESHOLDS_table.xlsx")

##Thresholds -----------------------------------------------------------------------------
##T2D ----------------------------------------------------------------------------------------
### 5%
UNITED_type2_all_genes_old_5PERC <- UNITED_type2_all_genes_old %>%
  summarise(totalover = sum(prob >= 0.05),
            ncasespickedup = sum(prob >= 0.05 & M ==1),
            PPV = (sum(prob >= 0.05 & M ==1)/sum(prob >= 0.05))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.05 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.05 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.05 & M == 0)/sum(prob < 0.05))*100,
            Sensitivity = (sum(prob >= 0.05 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.05 & M == 0)/sum(M==0))*100)
### 10%
UNITED_type2_all_genes_old_10PERC <- UNITED_type2_all_genes_old %>%
  summarise(totalover = sum(prob >= 0.1),
            ncasespickedup = sum(prob >= 0.1 & M ==1),
            PPV = (sum(prob >= 0.1 & M ==1)/sum(prob >= 0.1))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.1 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.1 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.1 & M == 0)/sum(prob < 0.1))*100,
            Sensitivity = (sum(prob >= 0.1 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.1 & M == 0)/sum(M==0))*100)

### 20%
UNITED_type2_all_genes_old_20PERC <- UNITED_type2_all_genes_old %>%
  summarise(totalover = sum(prob >= 0.2),
            ncasespickedup = sum(prob >= 0.2 & M ==1),
            PPV = (sum(prob >= 0.2 & M ==1)/sum(prob >= 0.2))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.2 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.2 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.2 & M == 0)/sum(prob < 0.2))*100,
            Sensitivity = (sum(prob >= 0.2 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.2 & M == 0)/sum(M==0))*100)

UNITED_type2_all_genes_old_THRESHOLDS <- rbind(UNITED_type2_all_genes_old_5PERC, UNITED_type2_all_genes_old_10PERC, UNITED_type2_all_genes_old_20PERC)
write_xlsx(UNITED_type2_all_genes_old_THRESHOLDS,"UNITED_type2_all_genes_old_THRESHOLDS_table.xlsx")
