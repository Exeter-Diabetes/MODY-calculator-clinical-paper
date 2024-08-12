#UNITED Paediatric PARTICIPANT CHARACTERISTICS TABLE --------------------------------------------

# load libraries
library(nimble)
library(rms)
library(tidyverse)
library(writexl)

## load functions --------------------------------------------
# needed for generating data 
source("data/create_data.R")
# load function needed for generating patient characteristics table
source("var_characteristics.R")

#load in UNITED type 1 datasets -----------------------------------------------------
## Just 3 genes ------------------------------------------------------------------------------
#M = MODY status (NA = not tested for MODY (i.e. T=1), 1 = tested for MODY (T=0) and positive, 0 = tested for MODY (T=0) and negative)
#C = c-peptide status (1 = UCPCR >= 0.2; 0 = UCPCR < 0.2)
#A = antibody status (1 = 1+ positive antibody, 0 = all antibodies tested negative)
#pardm = parent with diabetes (1 = 1+ parent affected, 0 = no parents affected)
#insoroha = currently treated with insulin or oha (tablets) (1= currently treated with either, 0 = not treated with either)
#T = biomarker status (1 = cpeptide negative (UCPCR < 0.2) or antibody positive (A =1), 0 = cpeptide positive (UCPCR >=0.2) AND antibody negative (A=0))
dataset.UNITED_type1p <- create_data(dataset = "united t1d pediatrics", commonmody = FALSE)

#need to change M=NA to M=0
dataset.UNITED_type1p <- dataset.UNITED_type1p %>%
  mutate(M = ifelse(is.na(M) == TRUE, 0, M))
#checked if worked: should have M=1 (n=7) & M=0 (n=1164)
table(dataset.UNITED_type1p$M)

#Produce characteristics tables ------------------------------------------------------------
## for early-insulin-treated ----------------------------------------------------------------
#create varlist (numeric variables of interest names)
varlist_T1D = c("agedx","agerec", "bmi", "hba1c")
#create varlist_cat (categorical variables of interest names)
varlist_cat_T1D = c("sex", "pardm", "insoroha", "C", "A", "T")
### for 3 genes ----------------------------------------------------------------------------------
#create table for UNITED T1D by MODY status - chosen to have numeric variables displayed as median [IQR]
var_characteristics(varlist = varlist_T1D, varlist_cat = varlist_cat_T1D, dataset = dataset.UNITED_type1p, numeric_option = "medianIQR", group = "M")
#save as
UNITED_T1Dp_table <- as.data.frame(summaryTable_GROUP_missing)
write_xlsx(UNITED_T1Dp_table,"UNITED_T1Dp_table.xlsx")



