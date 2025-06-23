##########################################################################################
#Table 1

#This code produces Table 1 and Supplementary Table 5 

#These tables contain the diagnostic & discriminatory criteria resulting 
#from using either a 5, 10, 15, 20, 25, or 30% testing threshold
#applied to the whole UNITED cohort used in this paper for both early- and not-
#early-insulin-treated participants

###############################################################################################
##Clinical utility table UNITED -----------------------------------------------------------------
#setwd("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper")
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
                                      id = TRUE, 
                                      biomarkers = "full") %>%
  ## if MODY testing missing, change to 0
  mutate(M = ifelse(is.na(M), 0, M),
         MODEL = "type1")
UNITED_type2_all_genes <- create_data(dataset = "united t2d", 
                                      commonmody = FALSE, 
                                      id = TRUE, 
                                      biomarkers = "full") %>%
  rename(C = UCPCRPosNegFinal, 
         A = AntibodyFINAL) %>%
  mutate(MODEL = "type2",
         A = ifelse(M == 1 & is.na(A), 0, A))


#load predictions ------------------------------------------------------------------------------------
#load T1D
predictions_UNITED_type1_with_T <- readRDS("Model_Predictions/predictions_dataset.UNITED_type1_all_genes_with_T.rds")
UNITED_type1 <- cbind(UNITED_type1_all_genes, predictions_UNITED_type1_with_T)
#load T2D
predictions_UNITED_type2 <- readRDS("Model_Predictions/predictions_dataset.UNITED_type2_all_genes_new.rds")
UNITED_type2 <- cbind(UNITED_type2_all_genes, predictions_UNITED_type2)


#Join datasets -----------------------------------------------------------------------
UNITED_joint <- full_join(UNITED_type1, UNITED_type2) %>%
  mutate(clinical = ifelse(MODEL == "type2" & pardm == 1 & agedx < 25, 1, 0))

#Load function to get diagnostic info at different thresholds ---------------------------------------
## Get threshold information -----------------------------------------------------------------
#Save as a table ---------------------------------------------------------------------
###UNITED pathway
UNITED_BIO <- UNITED_joint %>%
  summarise(totalover = sum((MODEL == "type1" & T == 0) | 
                              (MODEL == "type2"  & A == 0)),
            ncasespickedup = sum(((MODEL == "type1" & T == 0) | 
                                    (MODEL == "type2" & A == 0)) 
                                 & M ==1),
            PPV = (sum(((MODEL == "type1" & T == 0) | 
                          (MODEL == "type2"  & A == 0)) 
                       & M ==1)/sum((MODEL == "type1" & T == 0) | 
                                      (MODEL == "type2" 
                                        & A == 0)))*100,
            nmissedcases = sum(M==1) - sum(((MODEL == "type1" & T == 0) | 
                                             (MODEL == "type2" & A == 0)) 
                                           & M ==1),
            Missedcases = ((sum(M==1) - sum(((MODEL == "type1" & T == 0) | 
                                              (MODEL == "type2"  & A == 0)) 
                                            & M ==1))/sum(M==1))*100,
            NPV = (sum(((MODEL == "type1" & T == 1) | 
                          (MODEL == "type2" & A == 1)) & M == 0)/sum((MODEL == "type1" & T == 1) | 
                                                                       (MODEL == "type2" & A == 1)))*100,
            Sensitivity = (sum(((MODEL == "type1" & T == 0) | 
                                 (MODEL == "type2" & A == 0)) 
                               & M ==1)/sum(M==1))*100,
            Specificity = (sum(((MODEL == "type1" & T == 1) | 
                                  (MODEL == "type2" & A == 1)) & M == 0)/sum(M==0))*100)

### clinical 
UNITED_clinical <- UNITED_joint %>%
  summarise(totalover = sum(clinical == 1),
            ncasespickedup = sum(clinical == 1 & M ==1),
            PPV = (sum(clinical == 1 & M ==1)/sum(clinical == 1))*100,
            nmissedcases = sum(M==1) - sum(clinical == 1 & M ==1),
            Missedcases = ((sum(M==1) - sum(clinical == 1 & M ==1))/sum(M==1))*100,
            NPV = (sum(clinical == 0 & M == 0)/sum(clinical == 0))*100,
            Sensitivity = (sum(clinical == 1 & M ==1)/sum(M==1))*100,
            Specificity = (sum(clinical == 0 & M == 0)/sum(M==0))*100)
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


UNITED_compare_testing <- rbind(UNITED_5PERC,
                                UNITED_BIO,
                                UNITED_clinical)
write_xlsx(UNITED_compare_testing,"Supp_Table6.xlsx")



