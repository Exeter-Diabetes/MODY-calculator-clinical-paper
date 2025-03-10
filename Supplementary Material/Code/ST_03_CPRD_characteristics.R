##CPRD Characteristics Table

# load libraries
library(nimble)
library(rms)
library(tidyverse)
library(writexl)


#load functions -----------------------------------------------------------------------
source("var_characteristics.R")
#load data -----------------------------------------------------------------------------
source("pedro_mody_cohort.R")

#load in CPRD dataset
#M = MODY status (NA = not tested for MODY (i.e. T=1), 1 = tested for MODY (T=0) and positive, 0 = tested for MODY (T=0) and negative)
#C = c-peptide status (1 = UCPCR >= 0.2; 0 = UCPCR < 0.2)
#A = antibody status (1 = 1+ positive antibody, 0 = all antibodies tested negative)
#pardm = parent with diabetes (1 = 1+ parent affected, 0 = no parents affected)
#insoroha = currently treated with insulin or oha (tablets) (1= currently treated with either, 0 = not treated with either)
#T = biomarker status (1 = cpeptide negative (UCPCR < 0.2) or antibody positive (A =1), 0 = cpeptide positive (UCPCR >=0.2) AND antibody negative (A=0))
pedro_mody_cohort_local <- pedro_mody_cohort_local %>%
  mutate(pardm = ifelse(is.na(pardm), 0, pardm))

dataset.CPRD_type1 <- pedro_mody_cohort_local %>%
  filter(which_equation == "t1")
#load in UNITED type 2 dataset
dataset.CPRD_type2 <- pedro_mody_cohort_local %>%
  filter(which_equation == "t2")

dataset.CPRD_type1 <- dataset.CPRD_type1 %>%
  filter(ethnicity_5cat == "0")
#load in UNITED type 2 dataset
dataset.CPRD_type2 <- dataset.CPRD_type2 %>%
  filter(ethnicity_5cat == "0")

dataset.CPRD_type1 <- dataset.CPRD_type1 %>%
 drop_na(agedx, agerec, hba1c, bmi, pardm, sex)
#load in UNITED type 2 dataset
dataset.CPRD_type2 <- dataset.CPRD_type2 %>%
  drop_na(agedx, agerec, hba1c, bmi, pardm, sex, insoroha)

dataset.CPRD_type1_30 <- dataset.CPRD_type1 %>%
  filter(agedx <= 30)
dataset.CPRD_type2_30 <- dataset.CPRD_type2 %>%
  filter(agedx <= 30)

dataset.CPRD_type1_35 <- dataset.CPRD_type1 %>%
  filter(agedx <= 35)
dataset.CPRD_type2_35 <- dataset.CPRD_type2 %>%
  filter(agedx <= 35)





#PRODUCE TABLES --------------------------------------------------------------------------
#filter(!is.na(dm_diag_age) & !is.na(index_bmi_post_diag) & !is.na(index_hba1c_post_diag) & !is.na(gender) & !is.na(current_oha) & !is.na(fh_diabetes0) & !is.na(age_at_index))

#create varlist (numeric variables of interest names)
varlist = c("agedx","agerec", "bmi", "hba1c")
#create varlist_cat (categorical variables of interest names)
varlist_cat = c("sex", "pardm", "insoroha")

#sex = 1 is female, 0 is male


#create table for SA_MODY_5cat_full - chosen to have numeric variables displayed as meanCI
var_characteristics(varlist = varlist, varlist_cat = varlist_cat, dataset = dataset.CPRD_type1_35, numeric_option = "medianIQR")

#save as
CPRD_type1_table <- as.data.frame(summaryTable_count_num_miss)
write_xlsx(CPRD_type1_table,"Supplementary Material/Outputs/CPRD_type1_table.xlsx")

#create table for SA_MODY_5cat_full - chosen to have numeric variables displayed as meanCI
var_characteristics(varlist = varlist, varlist_cat = varlist_cat, dataset = dataset.CPRD_type2_35, numeric_option = "medianIQR")

#save as
CPRD_type2_table <- as.data.frame(summaryTable_count_num_miss)
write_xlsx(CPRD_type2_table,"Supplementary Material/Outputs/CPRD_type2_table.xlsx")