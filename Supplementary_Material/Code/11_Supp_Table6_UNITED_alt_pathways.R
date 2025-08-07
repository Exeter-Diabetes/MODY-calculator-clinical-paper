##########################################################################################
#Supplementary Table 6

#This code produces Supplementary Table 6 

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
##Apply Monogenice calculator predictions -----------------------------------------------
predictions_UNITED_type1_with_T <- readRDS("Model_Predictions/predictions_dataset.UNITED_type1_all_genes_with_T.rds") %>% 
  as.data.frame() %>%
  { rownames(.) <- NULL; . } %>%
  column_to_rownames(var = "id")
predictions_UNITED_type1_with_T <- predictions_UNITED_type1_with_T[as.character(c(UNITED_type1_all_genes$id)), ]
UNITED_type1 <- cbind(UNITED_type1_all_genes, predictions_UNITED_type1_with_T)
#load T2D
predictions_UNITED_type2 <- readRDS("Model_Predictions/predictions_dataset.UNITED_type2_all_genes_new.rds") %>% 
  as.data.frame() %>%
  { rownames(.) <- NULL; . } %>%
  column_to_rownames(var = "id")
predictions_UNITED_type2 <- predictions_UNITED_type2[as.character(c(UNITED_type2_all_genes$id)), ]
UNITED_type2 <- cbind(UNITED_type2_all_genes, predictions_UNITED_type2)

##Apply old calculator probabilities---------------------------------------------------
#from there use the separate regression equations = LogOR in separate equations
#less than 6 months
UNITED_type1 <- UNITED_type1 %>%
  mutate(
    old_logor = 1.8196 + (3.1404*pardm) - (0.0829*agerec) - (0.6598*hba1c) + 
      (0.1011*agedx) + (1.3131*sex), 
    old_prob = exp(old_logor)/(1 + exp(old_logor)),
    old_test = ifelse(old_prob >= 0.1, "Test", "Don't test"), 
    NHS_test = ifelse(T==0 & old_test == "Test", "Test", "Don't test")
  )

UNITED_type2 <- UNITED_type2 %>%
  mutate(
    old_logor = 19.28 - (0.3154*agedx) - (0.2324*bmi) - (0.6276*hba1c) + 
      (1.7473*pardm) - (0.0352*agerec) - (0.9952*insoroha) + (0.6943*sex), 
    old_prob = exp(old_logor)/(1 + exp(old_logor)),
    old_test = ifelse(old_prob >= 0.2, "Test", "Don't test"),
    NHS_test = ifelse(old_prob >= 0.2, "Test", "Don't test")
  )
#Join datasets -----------------------------------------------------------------------
UNITED_joint <- full_join(UNITED_type1, UNITED_type2) %>%
  mutate(clinical = ifelse(MODEL == "type2" & pardm == 1 & agedx < 25, 1, 0))

#Load function to get diagnostic info at different thresholds ---------------------------------------
## Get threshold information -----------------------------------------------------------------
#Save as a table ---------------------------------------------------------------------
###UNITED pathway
UNITED_BIO <- UNITED_joint %>%
  summarise(
    N = n(),
    nM = sum(M == 1),
    nNM = sum(M == 0),
    totalover = sum((MODEL == "type1" & T == 0) | (MODEL == "type2"  & A == 0)),
    totalunder = sum((MODEL == "type1" & T == 1) | (MODEL == "type2" & A == 1)),
    perc_testing = (totalover/N)*100,
    ncasespickedup = sum(((MODEL == "type1" & T == 0)| (MODEL == "type2" & A == 0)) & M ==1),
    ncontrolmissed = sum(((MODEL == "type1" & T == 1) | (MODEL == "type2" & A == 1)) & M == 0),
    PPV = (ncasespickedup/totalover)*100,
    nmissedcases = nM - ncasespickedup,
    Missedcases = (nmissedcases/nM)*100,
    NPV = (ncontrolmissed/totalunder)*100,
    Sens = (ncasespickedup/nM)*100,
    Spec = (ncontrolmissed/nNM)*100) 

### clinical 
UNITED_clinical <- UNITED_joint %>%
  summarise(
    N = n(),
    nM = sum(M == 1),
    nNM = sum(M == 0),
    totalover = sum(clinical == 1),
    totalunder = sum(clinical == 0),
    perc_testing = (totalover/N)*100,
    ncasespickedup = sum(clinical == 1 & M ==1),
    ncontrolmissed = sum(clinical == 0 & M ==0),
    PPV = (ncasespickedup/totalover)*100,
    nmissedcases = nM - ncasespickedup,
    Missedcases = (nmissedcases/nM)*100,
    NPV = (ncontrolmissed/totalunder)*100,
    Sens = (ncasespickedup/nM)*100,
    Spec = (ncontrolmissed/nNM)*100) 
### 5%
UNITED_5PERC <- UNITED_joint %>%
  summarise(
    N = n(),
    nM = sum(M == 1),
    nNM = sum(M == 0),
    totalover = sum(prob >= 0.05),
    totalunder = sum(prob < 0.05),
    perc_testing = (totalover/N)*100,
    ncasespickedup = sum(prob >= 0.05 & M ==1),
    ncontrolmissed = sum(prob < 0.05 & M == 0),
    PPV = (ncasespickedup/totalover)*100,
    nmissedcases = nM - ncasespickedup,
    Missedcases = (nmissedcases/nM)*100,
    NPV = (ncontrolmissed/totalunder)*100,
    Sens = (ncasespickedup/nM)*100,
    Spec = (ncontrolmissed/nNM)*100) 

### 5% MODY calculator
UNITED_5old <- UNITED_joint %>%
  summarise(
    N = n(),
    nM = sum(M == 1),
    nNM = sum(M == 0),
    totalover = sum(old_prob >= 0.05),
    totalunder = sum(old_prob < 0.05),
    perc_testing = (totalover/N)*100,
    ncasespickedup = sum(old_prob >= 0.05 & M ==1),
    ncontrolmissed = sum(old_prob < 0.05 & M == 0),
    PPV = (ncasespickedup/totalover)*100,
    nmissedcases = nM - ncasespickedup,
    Missedcases = (nmissedcases/nM)*100,
    NPV = (ncontrolmissed/totalunder)*100,
    Sens = (ncasespickedup/nM)*100,
    Spec = (ncontrolmissed/nNM)*100) 
#Using 10 and 20 % cut offs for MODY calculator
UNITED_old <- UNITED_joint %>%
  summarise(
    N = n(),
    nM = sum(M == 1),
    nNM = sum(M == 0),
    totalover = sum(old_test == "Test"),
    totalunder = sum(old_test == "Don't test"),
    perc_testing = (totalover/N)*100,
    ncasespickedup = sum(old_test == "Test" & M ==1),
    ncontrolmissed = sum(old_test == "Don't test" & M ==0),
    PPV = (ncasespickedup/totalover)*100,
    nmissedcases = nM - ncasespickedup,
    Missedcases = (nmissedcases/nM)*100,
    NPV = (ncontrolmissed/totalunder)*100,
    Sens = (ncasespickedup/nM)*100,
    Spec = (ncontrolmissed/nNM)*100) 
#Using 10 and 20 % cut offs for MODY calculator
#In C+, A- subset
UNITED_subset <- UNITED_joint %>%
  filter(T == 0) %>%
  summarise(
    N = n(),
    nM = sum(M == 1),
    nNM = sum(M == 0),
    totalover = sum(old_test == "Test"),
    totalunder = sum(old_test == "Don't test"),
    perc_testing = (totalover/N)*100,
    ncasespickedup = sum(old_test == "Test" & M ==1),
    ncontrolmissed = sum(old_test == "Don't test" & M ==0),
    PPV = (ncasespickedup/totalover)*100,
    nmissedcases = nM - ncasespickedup,
    Missedcases = (nmissedcases/nM)*100,
    NPV = (ncontrolmissed/totalunder)*100,
    Sens = (ncasespickedup/nM)*100,
    Spec = (ncontrolmissed/nNM)*100) 
#Using 5% cut offs for MODY calculator
#In C+, A- subset
UNITED_old_subset5 <- UNITED_joint %>%
  filter(T == 0) %>%
  summarise(
    N = n(),
    nM = sum(M == 1),
    nNM = sum(M == 0),
    totalover = sum(old_prob >= 0.05),
    totalunder = sum(old_prob < 0.05),
    perc_testing = (totalover/N)*100,
    ncasespickedup = sum(old_prob >= 0.05 & M ==1),
    ncontrolmissed = sum(old_prob < 0.05 & M == 0),
    PPV = (ncasespickedup/totalover)*100,
    nmissedcases = nM - ncasespickedup,
    Missedcases = (nmissedcases/nM)*100,
    NPV = (ncontrolmissed/totalunder)*100,
    Sens = (ncasespickedup/nM)*100,
    Spec = (ncontrolmissed/nNM)*100) 

UNITED_subset5 <- UNITED_joint %>%
  filter(T == 0) %>%
  summarise(
    N = n(),
    nM = sum(M == 1),
    nNM = sum(M == 0),
    totalover = sum(prob >= 0.05),
    totalunder = sum(prob < 0.05),
    perc_testing = (totalover/N)*100,
    ncasespickedup = sum(prob >= 0.05 & M ==1),
    ncontrolmissed = sum(prob < 0.05 & M == 0),
    PPV = (ncasespickedup/totalover)*100,
    nmissedcases = nM - ncasespickedup,
    Missedcases = (nmissedcases/nM)*100,
    NPV = (ncontrolmissed/totalunder)*100,
    Sens = (ncasespickedup/nM)*100,
    Spec = (ncontrolmissed/nNM)*100) 

#Using 10% cut offs for MODY calculator
#In early-insulin-treated individuals
UNITED_old_t1d <- UNITED_joint %>%
  filter(MODEL == "type1") %>%
  summarise(
    N = n(),
    nM = sum(M == 1),
    nNM = sum(M == 0),
    totalover = sum(old_prob >= 0.1),
    totalunder = sum(old_prob < 0.1),
    perc_testing = (totalover/N)*100,
    ncasespickedup = sum(old_prob >= 0.1 & M ==1),
    ncontrolmissed = sum(old_prob < 0.1 & M == 0),
    PPV = (ncasespickedup/totalover)*100,
    nmissedcases = nM - ncasespickedup,
    Missedcases = (nmissedcases/nM)*100,
    NPV = (ncontrolmissed/totalunder)*100,
    Sens = (ncasespickedup/nM)*100,
    Spec = (ncontrolmissed/nNM)*100) 

#Using 20% threshold for not-early-insulin-treated
#MODY calculator
UNITED_NHS_old_t2d <- UNITED_joint %>%
  filter(MODEL == "type2") %>%
  summarise(
    N = n(),
    nM = sum(M == 1),
    nNM = sum(M == 0),
    totalover = sum(old_prob >= 0.2),
    totalunder = sum(old_prob < 0.2),
    perc_testing = (totalover/N)*100,
    ncasespickedup = sum(old_prob >= 0.2 & M ==1),
    ncontrolmissed = sum(old_prob < 0.2 & M == 0),
    PPV = (ncasespickedup/totalover)*100,
    nmissedcases = nM - ncasespickedup,
    Missedcases = (nmissedcases/nM)*100,
    NPV = (ncontrolmissed/totalunder)*100,
    Sens = (ncasespickedup/nM)*100,
    Spec = (ncontrolmissed/nNM)*100) 

#Using 25% threshold across all MODY calculator
UNITED_25old <- UNITED_joint %>%
  summarise(
    N = n(),
    nM = sum(M == 1),
    nNM = sum(M == 0),
    totalover = sum(old_prob >= 0.25),
    totalunder = sum(old_prob < 0.25),
    perc_testing = (totalover/N)*100,
    ncasespickedup = sum(old_prob >= 0.25 & M ==1),
    ncontrolmissed = sum(old_prob < 0.25 & M == 0),
    PPV = (ncasespickedup/totalover)*100,
    nmissedcases = nM - ncasespickedup,
    Missedcases = (nmissedcases/nM)*100,
    NPV = (ncontrolmissed/totalunder)*100,
    Sens = (ncasespickedup/nM)*100,
    Spec = (ncontrolmissed/nNM)*100) 

UNITED_NHS_old_t1d <- UNITED_joint %>%
  filter(MODEL == "type1") %>%
  summarise(
    N = n(),
    nM = sum(M == 1),
    nNM = sum(M == 0),
    totalover = sum(NHS_test == "Test"),
    totalunder = sum(NHS_test == "Don't test"),
    perc_testing = (totalover/N)*100,
    ncasespickedup = sum(NHS_test == "Test" & M ==1),
    ncontrolmissed = sum(NHS_test == "Don't test" & M ==0),
    PPV = (ncasespickedup/totalover)*100,
    nmissedcases = nM - ncasespickedup,
    Missedcases = (nmissedcases/nM)*100,
    NPV = (ncontrolmissed/totalunder)*100,
    Sens = (ncasespickedup/nM)*100,
    Spec = (ncontrolmissed/nNM)*100) 

UNITED_NHS_old <- UNITED_joint %>%
  summarise(
    N = n(),
    nM = sum(M == 1),
    nNM = sum(M == 0),
    totalover = sum(NHS_test == "Test"),
    totalunder = sum(NHS_test == "Don't test"),
    perc_testing = (totalover/N)*100,
    ncasespickedup = sum(NHS_test == "Test" & M ==1),
    ncontrolmissed = sum(NHS_test == "Don't test" & M ==0),
    PPV = (ncasespickedup/totalover)*100,
    nmissedcases = nM - ncasespickedup,
    Missedcases = (nmissedcases/nM)*100,
    NPV = (ncontrolmissed/totalunder)*100,
    Sens = (ncasespickedup/nM)*100,
    Spec = (ncontrolmissed/nNM)*100) 

#Merge into one table --------------------------------------------------------------------
UNITED_compare_testing <- rbind(
  UNITED_5PERC,
  UNITED_BIO,
  UNITED_clinical, 
  UNITED_5old, 
  #UNITED_old,
  UNITED_NHS_old,
  #UNITED_NHS_old_t1d,
  #UNITED_NHS_old_t2d,
  #UNITED_25old,
  #UNITED_old_t1d,
  UNITED_subset5,
  UNITED_old_subset5#,
  #UNITED_subset
  ) %>%
  mutate(
    `% individuals tested at threshold` = paste0(
      round(perc_testing,1), " (", totalover, "/", N, ")"),
    `% monogenic diabetes cases missed` = paste0(
      round(Missedcases,1), " (", nmissedcases, "/", nM, ")"),
    `% positive test rate (PPV)` = paste0(
      round(PPV,1), " (",ncasespickedup, "/", totalover, ")"),
    `% negative test rate (NPV)` = paste0(
      round(NPV,1), " (", ncontrolmissed, "/", totalunder, ")"),
    Sensitivity = paste0(
      round(Sens,1), " (", ncasespickedup, "/", nM,")"), 
    Specificity = paste0(
      round(Spec,1), " (", ncontrolmissed, "/", nNM, ")")
  ) %>%
  select(
    `% individuals tested at threshold`,
    `% monogenic diabetes cases missed`,
    `% positive test rate (PPV)`,
    `% negative test rate (NPV)`,
    Sensitivity, 
    Specificity
  )
write_xlsx(UNITED_compare_testing,"Supp_Table6.xlsx")



