#####################################################################################
#Supplementary Table 4 - Characteristics table of external validation dataset
##########################################################################################
#SET WORKING DIRECTORY -----------------------------------------------------------------------------------
setwd("~/PhD/CLINICAL MODY/MODY-calculator-clinical-paper")

#load libraries
library(tidyverse)
library(writexl)
library(pROC)


#load functions
source("Functions/var_characteristics_1.R")
source("Data/create_data.R")

#load data
# MyDiabetes ------------------------------------------------------------------------------------
#load T1D
load("Data/MY_T1D.RData")
#load T2D 
load("Data/MY_T2D.RData")
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
UNITED_p <- create_data(dataset = "united t1d pediatrics", 
                        commonmody = FALSE, 
                        id = TRUE)

#need to change M=NA to M=0
UNITED_p <- UNITED_p %>%
  mutate(M = ifelse(is.na(M) == TRUE, 0, M)) %>%
  filter(!is.na(T))
#checked if worked: should have M=1 (n=7) & M=0 (n=1164)
table(UNITED_p$M)
table(UNITED_p$M, UNITED_p$Gene)

#T2D
UNITED_type2p <- UNITED_p %>%
  filter(tti == "" | tti == "Greater than 12 months")
table(UNITED_type2p$M)

#T1D
UNITED_type1p <- UNITED_p %>%
  filter(!(id %in% UNITED_type2p$id))
table(UNITED_type1p$M)

# External datasets -----------------------------------------------------------------------
## Early-insulin-treated
MYDIABETES_type1 <- MY_T1D %>%
  mutate(study = "MYDIABETES",
         durationfinal = agerec - agedx) %>%
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
         biomark_status,
         durationfinal) 

UNITED_type1p <- UNITED_type1p %>%
  mutate(study = "UNITED paediatric") %>%
  filter(!is.na(T))
dataset_type1 <- full_join(MYDIABETES_type1, 
                           UNITED_type1p, 
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
                                  "Gene", 
                                  "durationfinal")) %>%
  mutate(M = ifelse(is.na(M), 0, M), 
         Gene = ifelse(Gene == "", NA, Gene), 
         biomark_status = ifelse(is.na(biomark_status), 
                                 ifelse(C == 1 & A == 0, 0, 1),
                                 biomark_status),
         type = ifelse(M == 1,
                       "Mody",
                       "T1D")) 
# Not-early-insulin-treated
MYDIABETES_type2 <- MY_T2D %>%
  mutate(study = "MYDIABETES",
         durationfinal = agerec - agedx) %>%
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
         biomark_status, 
         durationfinal) 

UNITED_type2p <- UNITED_type2p %>%
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
         Gene,
         durationfinal) 

dataset_type2 <- full_join(MYDIABETES_type2, 
                           UNITED_type2p, 
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
                                  "Gene",
                                  "durationfinal")) %>%
  mutate(M = ifelse(is.na(M), 0, M), 
         Gene = ifelse(Gene == "", NA, Gene), 
         biomark_status = ifelse(is.na(biomark_status), 
                                 ifelse(C == 1 & A == 0, 0, 1),
                                 biomark_status),
         type = ifelse(M == 1,
                       "Mody",
                       "T2D")) 

external_data <- full_join(dataset_type1, dataset_type2)
mydiabetes <- full_join(MYDIABETES_type1, MYDIABETES_type2)

###########################################################################################
#1. CHARACTERISTICS TABLES
#create varlist (numeric variables of interest names)
varlist = c("agedx",
            "agerec", 
            "bmi", 
            "hba1c",
            "durationfinal")
#create varlist_cat (categorical variables of interest names)
varlist_cat = c("Gene", 
                "sex", 
                "pardm", 
                "insoroha", 
                "C", 
                "A", 
                "biomark_status")

# Early-insulin-treated ------------------------------------------------------------------------------------------
#create table for external data - early-insulin-treated by MODY status   
#- chosen to have numeric variables displayed as median [IQR]
var_characteristics(varlist = varlist, 
                    varlist_cat = varlist_cat, 
                    dataset = dataset_type1, 
                    numeric_option = "medianIQR", 
                    group = "M",
                    p_value_testing = FALSE,
                    table_name = "Supplementary_Material/Outputs/Supp_Table4_T1D_table")


# Not-early-insulin-treated  -----------------------------------------------------------------------------------------
#create table for external data - not-early-insulin-treated by MODY status 
#- chosen to have numeric variables displayed as median [IQR]
var_characteristics(varlist = varlist, 
                    varlist_cat = varlist_cat, 
                    dataset = dataset_type2, 
                    numeric_option = "medianIQR", 
                    group = "M",
                    p_value_testing = FALSE,
                    table_name = "Supplementary_Material/Outputs/Supp_Table4_T2D_table")


#Joint ---------------------------------------------------------------------------------
#create table for UNITED T1D by MODY status - chosen to have numeric variables displayed as median [IQR]
var_characteristics(varlist = varlist, 
                    varlist_cat = varlist_cat, 
                    dataset = external_data, 
                    numeric_option = "medianIQR", 
                    group = "type",
                    p_value_testing = FALSE,
                    table_name = "Supplementary_Material/Outputs/Supp_Table4_table")

external_data %>%
  summarise(
    median = median(durationfinal, na.rm = TRUE),
    Q1 = quantile(durationfinal, probs = 0.25, na.rm = TRUE),
    Q3 = quantile(durationfinal, probs = 0.75, na.rm = TRUE)
  )

table(external_data$type, external_data$Gene, useNA = "ifany")

UNITED_p %>%
  summarise(
    median = round(median(durationfinal, na.rm = TRUE),2),
    Q1 = round(quantile(durationfinal, probs = 0.25, na.rm = TRUE),2),
    Q3 = round(quantile(durationfinal, probs = 0.75, na.rm = TRUE),2)
  )
mydiabetes %>%
  summarise(
    median = round(median(durationfinal, na.rm = TRUE),2),
    Q1 = round(quantile(durationfinal, probs = 0.25, na.rm = TRUE),2),
    Q3 = round(quantile(durationfinal, probs = 0.75, na.rm = TRUE),2)
  )
