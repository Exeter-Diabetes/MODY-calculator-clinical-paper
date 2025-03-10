#Case-control characteristics table

# load libraries
library(nimble)
library(rms)
library(tidyverse)
library(writexl)

# load functions needed for generating data
source("data/create_data.R")

#load in Case-control type 1 dataset
#mody = MODY status (NA = not tested for MODY (i.e. T=1), 1 = tested for MODY (T=0) and positive, 0 = tested for MODY (T=0) and negative)
#pardm = parent with diabetes (1 = 1+ parent affected, 0 = no parents affected)
#insoroha = currently treated with insulin or oha (tablets) (1= currently treated with either, 0 = not treated with either)
dataset.cc_type1 <- create_data(dataset = "case-control t1d")
#load in UNITED type 2 dataset
dataset.cc_type2 <- create_data(dataset = "case-control t2d")

# load function needed for generating patient characteristics table
source("var_characteristics.R")

#create varlist (numeric variables of interest names)
varlist = c("agedx","agerec", "bmi", "hba1c")
#create varlist_cat (categorical variables of interest names)
varlist_cat_T1D = c("sex", "pardm")
#create table for UNITED T1D by MODY status - chosen to have numeric variables displayed as median [IQR]
var_characteristics(varlist = varlist, varlist_cat = varlist_cat_T1D, dataset = dataset.cc_type1, numeric_option = "medianIQR", group = "mody")

#save as
CC_T1D_table <- as.data.frame(summaryTable_GROUP_missing)


#create varlist_cat (categorical variables of interest names)
varlist_cat_T2D = c("sex", "pardm", "insoroha")
#create table for UNITED T1D by MODY status - chosen to have numeric variables displayed as median [IQR]
var_characteristics(varlist = varlist, varlist_cat = varlist_cat_T2D, dataset = dataset.cc_type2, numeric_option = "medianIQR", group = "mody")

#save as
CC_T2D_table <- as.data.frame(summaryTable_GROUP_missing)

write_xlsx(CC_T1D_table,"Supplementary Material/Outputs/CC_T1D_table.xlsx")
write_xlsx(CC_T2D_table,"Supplementary Material/Outputs/CC_T2D_table.xlsx")
