#UNITED PARTICIPANT CHARACTERISTICS TABLE

# load libraries
library(nimble)
library(rms)
library(tidyverse)

## load functions needed for generating data ----
source("data/create_data.R")

#load in UNITED type 1 dataset
#M = MODY status (NA = not tested for MODY (i.e. T=1), 1 = tested for MODY (T=0) and positive, 0 = tested for MODY (T=0) and negative)
#C = c-peptide status (1 = UCPCR >= 0.2; 0 = UCPCR < 0.2)
#A = antibody status (1 = 1+ positive antibody, 0 = all antibodies tested negative)
#pardm = parent with diabetes (1 = 1+ parent affected, 0 = no parents affected)
#insoroha = currently treated with insulin or oha (tablets) (1= currently treated with either, 0 = not treated with either)
#T = biomarker status (1 = cpeptide negative (UCPCR < 0.2) or antibody positive (A =1), 0 = cpeptide positive (UCPCR >=0.2) AND antibody negative (A=0))
dataset.UNITED_type1 <- create_data(dataset = "united t1d", biomarkers = "full")
#load in UNITED type 2 dataset
dataset.UNITED_type2 <- create_data(dataset = "united t2d")

#need to change M=NA to M=0
dataset.UNITED_type1 <- dataset.UNITED_type1 %>%
  mutate(M = ifelse(is.na(M) == TRUE, 0, M))
#checked if worked: should have M=1 (n=7) & M=0 (n=1164)
table(dataset.UNITED_type1$M)

# load function needed for generating patient characteristics table
source("var_characteristics.R")

#create varlist (numeric variables of interest names)
varlist_T1D = c("agedx","agerec", "bmi", "hba1c")
#create varlist_cat (categorical variables of interest names)
varlist_cat_T1D = c("sex", "pardm", "insoroha", "C", "A", "T")
#create table for UNITED T1D by MODY status - chosen to have numeric variables displayed as median [IQR]
var_characteristics(varlist = varlist_T1D, varlist_cat = varlist_cat_T1D, dataset = dataset.UNITED_type1, numeric_option = "medianIQR", group = "M")

#save as
UNITED_T1D_table <- as.data.frame(summaryTable_GROUP_missing)

#create varlist (numeric variables of interest names)
varlist_T2D = c("agedx","agerec", "bmi", "hba1c")
#create varlist_cat (categorical variables of interest names)
varlist_cat_T2D = c("sex", "pardm", "insoroha")
#create table for UNITED T1D by MODY status - chosen to have numeric variables displayed as median [IQR]
var_characteristics(varlist = varlist_T2D, varlist_cat = varlist_cat_T2D, dataset = dataset.UNITED_type2, numeric_option = "medianIQR", group = "M")

#save as
UNITED_T2D_table <- as.data.frame(summaryTable_GROUP_missing)

library(writexl)
write_xlsx(UNITED_T1D_table,"UNITED_T1D_table.xlsx")
write_xlsx(UNITED_T2D_table,"UNITED_T2D_table.xlsx")

#Finding the mean MODY prob for M+ AND M- participants
#load predictions
predictions_dataset.UNITED_type2_new <- readRDS("~/model_predictions/predictions_dataset.UNITED_type2_new.rds")
#bind them
UNITED_type2 <- cbind(dataset.UNITED_type2,predictions_dataset.UNITED_type2_new)

UNITED_type2 %>%
  filter(M==0) %>%
  summarise(mean = mean(prob))

UNITED_type2 %>%
  filter(M==1) %>%
  summarise(mean = mean(prob))

