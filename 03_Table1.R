##########################################################################################
#Table 1

#This code produces Table 1 and Supplementary Table 5 

#These tables contain the diagnostic & discriminatory criteria resulting 
#from using either a 5, 10, 15, 20, 25, or 30% testing threshold
#applied to the whole UNITED cohort used in this paper for both early- and not-
#early-insulin-treated participants

###############################################################################################
##Clinical utility table UNITED -----------------------------------------------------------------
setwd("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper")
#load libraries --------------------------------------------------------------------------------
library(tidyverse)
library(pROC)
library(PRROC)
library(writexl)

#load functions
source("Data/create_data.R")
## Extended table --------------------------------------------------------------------------
## All genes ---------------------------------------------------------------------------------------
UNITED_type1_all_genes <- create_data(dataset = "united t1d", 
                                      commonmody = FALSE, 
                                      id = TRUE) %>%
  ## if MODY testing missing, change to 0
  mutate(M = ifelse(is.na(M), 0, M))
UNITED_type2_all_genes <- create_data(dataset = "united t2d", 
                                      commonmody = FALSE, 
                                      id = TRUE)


#load predictions ------------------------------------------------------------------------------------
#load T1D
predictions_UNITED_type1_with_T <- readRDS("Model_Predictions/predictions_dataset.UNITED_type1_all_genes_with_T.rds")
UNITED_type1 <- cbind(dataset.UNITED_type1_all_genes, predictions_UNITED_type1_with_T)
#load T2D
predictions_UNITED_type2 <- readRDS("Model_Predictions/predictions_dataset.UNITED_type2_all_genes_new.rds")
UNITED_type2 <- cbind(dataset.UNITED_type2_all_genes, predictions_UNITED_type2)


#Join datasets -----------------------------------------------------------------------
UNITED_joint <- full_join(UNITED_type1, UNITED_type2)

#Load function to get diagnostic info at different thresholds ---------------------------------------
calculate_thresholds_diagnostics <- function(response, prediction, unique = FALSE) {
  
  ### threshold number to test
  thresholds <- seq(0, 1, 0.0001)
  
  matrix_thresholds <- matrix(0, nrow = length(thresholds), ncol = 8)
  matrix_thresholds[,1] <- thresholds
  
  preds1 <- data.frame(
    prediction = prediction,
    response = response
  )
  
  for (i in 1:length(thresholds)) {
    # original dataset
    original_data <- preds1
    
    # checking which patients are above the threshold
    interim_above <- preds1 %>%
      mutate(threshold_up = ifelse(prediction > thresholds[i], 1, 0)) %>%
      filter(threshold_up == 1)
    
    interim_below <- preds1 %>%
      mutate(threshold_up = ifelse(prediction > thresholds[i], 1, 0)) %>%
      filter(threshold_up == 0)
    
    # calculating the sensitivity of the MODY cases
    matrix_thresholds[i, 2] <- length(which(interim_above$response == 1))/length(which(original_data$response == 1)) * 100
    
    # calculating number to test
    if (length(which(interim_above$response == 1)) > 0) {
      matrix_thresholds[i, 3] <- ceiling(nrow(interim_above)/length(which(interim_above$response == 1)))
    } else {matrix_thresholds[i, 3] <- as.numeric(NA)}
    
    # Pick-up rate
    matrix_thresholds[i, 4] <- length(which(interim_above$response == 1))/nrow(interim_above) * 100
    
    # calculating specificity
    matrix_thresholds[i, 5] <- length(which(interim_below$response == 0))/length(which(original_data$response == 0)) * 100
    
    # Perc of patients tested
    matrix_thresholds[i, 6] <- (nrow(interim_above)/nrow(original_data))*100
    
    # number of MODY cases missed
    matrix_thresholds[i, 7] <- interim_below %>%
      select(response) %>%
      sum(na.rm = TRUE) %>%
      unlist()
    
    # number of patients tested
    matrix_thresholds[i, 8] <- interim_above %>%
      nrow() %>%
      unlist()
    
  }
  
  ## making the matrix a data.frame()
  matrix_thresholds <- matrix_thresholds %>% 
    as.data.frame() %>%
    set_names(c("Thresholds", "Sensitivity", "NTT", "Pick-up rate", "Specificity", "N-Tested", "Cases Missed", "Patients Tested")) %>%
    arrange(desc(Thresholds))
  
  if (unique == TRUE) {
    ## select unique combinations of sensitivity and NTT (only the first occurance)
    matrix_thresholds <- matrix_thresholds %>%
      slice(which(duplicated(matrix_thresholds %>% select(-Thresholds, -`Pick-up rate`)) == FALSE))
  }
  
  return(matrix_thresholds)
  
}


## Get threshold information -----------------------------------------------------------------
thresholds_UNITED <- calculate_thresholds_diagnostics(UNITED_joint$M, UNITED_joint$prob)

### 5%
thresholds_UNITED %>%
  filter(Thresholds == 0.05) %>% arrange(Thresholds) %>% head()


### 10%
thresholds_UNITED %>%
  filter(Thresholds == 0.1) %>% arrange(Thresholds)  %>% head()

### 15%
thresholds_UNITED %>%
  filter(Thresholds == 0.15) %>% arrange(Thresholds)  %>% head()


### 20%
thresholds_UNITED %>%
  filter(Thresholds == 0.2) %>% arrange(Thresholds)  %>% head()



### 25%
thresholds_UNITED %>%
  filter(Thresholds == 0.25) %>% arrange(Thresholds)  %>% head()



### 30%
thresholds_UNITED %>%
  filter(Thresholds == 0.3) %>% arrange(Thresholds)  %>% head()

#Save as a table ---------------------------------------------------------------------
### 5%
UNITED_5PERC <- UNITED_joint %>%
  summarise(totalover = sum(prob >= 0.05),
            ncasespickedup = sum(prob >= 0.05 & M ==1),
            PPV = (sum(prob >= 0.05 & M ==1)/sum(prob >= 0.05))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.05 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.05 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.05 & M == 0)/sum(prob < 0.05))*100,
            Sensitivity = (sum(prob >= 0.05 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.05 & M == 0)/sum(M==0))*100)
### 10%
UNITED_10PERC <- UNITED_joint %>%
  summarise(totalover = sum(prob >= 0.1),
            ncasespickedup = sum(prob >= 0.1 & M ==1),
            PPV = (sum(prob >= 0.1 & M ==1)/sum(prob >= 0.1))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.1 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.1 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.1 & M == 0)/sum(prob < 0.1))*100,
            Sensitivity = (sum(prob >= 0.1 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.1 & M == 0)/sum(M==0))*100)

### 15%
UNITED_15PERC <- UNITED_joint %>%
  summarise(totalover = sum(prob >= 0.15),
            ncasespickedup = sum(prob >= 0.15 & M ==1),
            PPV = (sum(prob >= 0.15 & M ==1)/sum(prob >= 0.15))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.15 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.15 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.15 & M == 0)/sum(prob < 0.15))*100,
            Sensitivity = (sum(prob >= 0.15 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.15 & M == 0)/sum(M==0))*100)

### 20%
UNITED_20PERC <- UNITED_joint %>%
  summarise(totalover = sum(prob >= 0.2),
            ncasespickedup = sum(prob >= 0.2 & M ==1),
            PPV = (sum(prob >= 0.2 & M ==1)/sum(prob >= 0.2))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.2 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.2 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.2 & M == 0)/sum(prob < 0.2))*100,
            Sensitivity = (sum(prob >= 0.2 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.2 & M == 0)/sum(M==0))*100)


### 25%
UNITED_25PERC <- UNITED_joint %>%
  summarise(totalover = sum(prob >= 0.25),
            ncasespickedup = sum(prob >= 0.25 & M ==1),
            PPV = (sum(prob >= 0.25 & M ==1)/sum(prob >= 0.25))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.25 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.25 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.25 & M == 0)/sum(prob < 0.25))*100,
            Sensitivity = (sum(prob >= 0.25 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.25 & M == 0)/sum(M==0))*100)

### 30%
UNITED_30PERC <- UNITED_joint %>%
  summarise(totalover = sum(prob >= 0.3),
            ncasespickedup = sum(prob >= 0.3 & M ==1),
            PPV = (sum(prob >= 0.3 & M ==1)/sum(prob >= 0.3))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.3 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.3 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.3 & M == 0)/sum(prob < 0.3))*100,
            Sensitivity = (sum(prob >= 0.3 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.3 & M == 0)/sum(M==0))*100)



UNITED_THRESHOLDS_Table1 <- rbind(UNITED_5PERC,
                                  UNITED_10PERC,
                                  UNITED_20PERC,
                                  UNITED_30PERC)
write_xlsx(UNITED_THRESHOLDS_Table1,"Table1.xlsx")



