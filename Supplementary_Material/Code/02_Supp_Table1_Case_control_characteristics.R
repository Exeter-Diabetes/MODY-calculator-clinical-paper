#Case-control characteristics table

# load libraries
library(tidyverse)
library(writexl)

# load functions needed for generating data
source("Data/create_data.R")
# load function needed for generating patient characteristics table
source("Functions/var_characteristics_1.R")

#load in Case-control type 1 dataset
#mody = MODY status (NA = not tested for MODY (i.e. T=1), 
#                     1 = tested for MODY (T=0) and positive, 
#                     0 = tested for MODY (T=0) and negative)
#pardm = parent with diabetes (1 = 1+ parent affected, 
#                              0 = no parents affected)
#insoroha = currently treated with insulin or oha (tablets) 
#           (1= currently treated with either, 
#            0 = not treated with either)

cc_type1 <- create_data(dataset = "case-control t1d")
#load in UNITED type 2 dataset
cc_type2 <- create_data(dataset = "case-control t2d")


#create varlist (numeric variables of interest names)
varlist = c("agedx",
            "agerec", 
            "bmi", 
            "hba1c")
#create varlist_cat (categorical variables of interest names)
varlist_cat_T1D = c("sex", 
                    "pardm")
#create table for UNITED T1D by MODY status - chosen to have numeric variables displayed as median [IQR]
var_characteristics(varlist = varlist, 
                    varlist_cat = varlist_cat_T1D, 
                    dataset = cc_type1, 
                    numeric_option = "medianIQR", 
                    group = "mody",
                    p_value_testing = FALSE,
                    table_name = "Supplementary_Material/Outputs/CC_T1D_table")


#create varlist_cat (categorical variables of interest names)
varlist_cat_T2D = c("sex", 
                    "pardm", 
                    "insoroha")
#create table for UNITED T1D by MODY status - chosen to have numeric variables displayed as median [IQR]
var_characteristics(varlist = varlist, 
                    varlist_cat = varlist_cat_T2D, 
                    dataset = cc_type2, 
                    numeric_option = "medianIQR", 
                    group = "mody",
                    p_value_testing = FALSE,
                    table_name = "Supplementary_Material/Outputs/CC_T2D_table")


