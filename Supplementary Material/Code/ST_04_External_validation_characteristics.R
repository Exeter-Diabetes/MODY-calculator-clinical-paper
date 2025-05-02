#####################################################################################
#Supplementary Table 4 - Characteristics table of external validation dataset
##########################################################################################
#SET WORKING DIRECTORY -----------------------------------------------------------------------------------
setwd("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper")

#load libraries
library(tidyverse)
library(writexl)
library(pROC)
#library(ggbeeswarm)

#load functions
source("var_characteristics.R")
source("data/create_data.R")

#load data
# MyDiabetes ------------------------------------------------------------------------------------
#load T1D
load("data/MY_T1D.RData")
#load T2D 
load("data/MY_T2D.RData")
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
dataset.UNITED_p <- create_data(dataset = "united t1d pediatrics", 
                                commonmody = FALSE, 
                                id = TRUE)

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
  select(MY_ID, 
         study, 
         agerec, 
         agedx, 
         sex, 
         bmi, 
         pardm, 
         insoroha, 
         hba1c, 
         C, 
         A, 
         M, 
         Gene, 
         biomark_status) 

dataset.UNITED_type1p <- dataset.UNITED_type1p %>%
  mutate(study = "UNITED paediatric") %>%
  filter(!is.na(T))
dataset_type1 <- full_join(dataset.MYDIABETES_type1, 
                           dataset.UNITED_type1p, 
                           by = c("study",
                                  "agerec", 
                                  "agedx", 
                                  "sex", 
                                  "bmi", 
                                  "pardm", "
                                  insoroha", 
                                  "hba1c", 
                                  "C", 
                                  "A", 
                                  "M", 
                                  "Gene")) %>%
  mutate(M = ifelse(is.na(M), 0, M), 
         Gene = ifelse(Gene == "", NA, Gene), 
         biomark_status = ifelse(is.na(biomark_status), 
                                 ifelse(C == 1 & A == 0, 0, 1),
                                 biomark_status),
         type = ifelse(M == 1,
                       "Mody",
                       "T1D")) 
# Not-early-insulin-treated
dataset.MYDIABETES_type2 <- MY_T2D %>%
  mutate(study = "MYDIABETES") %>%
  select(MY_ID, 
         study, 
         agerec, 
         agedx, 
         sex, 
         bmi, 
         pardm, 
         insoroha, 
         hba1c, 
         C, 
         A, 
         M, 
         Gene, 
         biomark_status) 

dataset.UNITED_type2p <- dataset.UNITED_type2p %>%
  mutate(study = "UNITED paediatric") %>%
  select(id, 
         study, 
         agerec, 
         agedx, 
         sex, 
         bmi, 
         pardm, 
         insoroha, 
         hba1c, 
         C, 
         A, 
         M, 
         Gene) 

dataset_type2 <- full_join(dataset.MYDIABETES_type2, 
                           dataset.UNITED_type2p, 
                           by = c("study",
                                  "agerec", 
                                  "agedx", 
                                  "sex", 
                                  "bmi", 
                                  "pardm", 
                                  "insoroha", 
                                  "hba1c", 
                                  "C", 
                                  "A", 
                                  "M", 
                                  "Gene")) %>%
  mutate(M = ifelse(is.na(M), 0, M), 
         Gene = ifelse(Gene == "", NA, Gene), 
         biomark_status = ifelse(is.na(biomark_status), 
                                 ifelse(C == 1 & A == 0, 0, 1),
                                 biomark_status),
         type = ifelse(M == 1,
                       "Mody",
                       "T2D")) 

external_data <- full_join(dataset_type1, dataset_type2)

###########################################################################################
#1. CHARACTERISTICS TABLES
#create varlist (numeric variables of interest names)
varlist = c("agedx","agerec", "bmi", "hba1c")
#create varlist_cat (categorical variables of interest names)
varlist_cat = c("Gene", "sex", "pardm", "insoroha", "C", "A", "biomark_status")

# Early-insulin-treated ------------------------------------------------------------------------------------------
#create table for external data - early-insulin-treated by MODY status   
#- chosen to have numeric variables displayed as median [IQR]
var_characteristics(varlist = varlist, 
                    varlist_cat = varlist_cat, 
                    dataset = dataset_type1, 
                    numeric_option = "medianIQR", 
                    group = "M")
#save as
externalval_t1d_character_table <- as.data.frame(summaryTable_GROUP_missing)
write_xlsx(externalval_t1d_character_table,
           "Supplementary Material/Outputs/ST_04_EARLY_INSULIN_TREATED_table.xlsx")

# Not-early-insulin-treated  -----------------------------------------------------------------------------------------
#create table for external data - not-early-insulin-treated by MODY status 
#- chosen to have numeric variables displayed as median [IQR]
var_characteristics(varlist = varlist, 
                    varlist_cat = varlist_cat, 
                    dataset = dataset_type2, 
                    numeric_option = "medianIQR", 
                    group = "M")
#save as
externalval_t2d_character_table <- as.data.frame(summaryTable_GROUP_missing)
write_xlsx(externalval_t2d_character_table,
           "Supplementary Material/Outputs/ST_04_NOT_EARLY_INSULIN_TREATED_table.xlsx")

#Joint ---------------------------------------------------------------------------------
#create table for UNITED T1D by MODY status - chosen to have numeric variables displayed as median [IQR]
var_characteristics(varlist = varlist, 
                    varlist_cat = varlist_cat, 
                    dataset = external_data, 
                    numeric_option = "medianIQR", 
                    group = "type")
#save as
externalval_character_table <- as.data.frame(summaryTable_GROUP_missing)
write_xlsx(externalval_character_table,"Supplementary Material/Outputs/ST_04_table.xlsx")
