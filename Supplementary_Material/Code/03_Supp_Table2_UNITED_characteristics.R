#UNITED PARTICIPANT CHARACTERISTICS TABLE --------------------------------------------

# load libraries
library(tidyverse)
library(writexl)

## load functions --------------------------------------------
# needed for generating data 
source("Data/create_data.R")
# load function needed for generating patient characteristics table
source("Functions/var_characteristics_1.R")

#load in UNITED type 1 datasets -----------------------------------------------------
## Just 3 genes ------------------------------------------------------------------------------
#M = MODY status (NA = not tested for MODY (i.e. T=1), 1 = tested for MODY (T=0) and positive, 0 = tested for MODY (T=0) and negative)
#C = c-peptide status (1 = UCPCR >= 0.2; 0 = UCPCR < 0.2)
#A = antibody status (1 = 1+ positive antibody, 0 = all antibodies tested negative)
#pardm = parent with diabetes (1 = 1+ parent affected, 0 = no parents affected)
#insoroha = currently treated with insulin or oha (tablets) (1= currently treated with either, 0 = not treated with either)
#T = biomarker status (1 = cpeptide negative (UCPCR < 0.2) or antibody positive (A =1), 0 = cpeptide positive (UCPCR >=0.2) AND antibody negative (A=0))
UNITED_type1 <- create_data(dataset = "united t1d",
                            commonmody = FALSE) %>%
  ## if MODY testing missing, change to 0
  mutate(M = ifelse(is.na(M), 0, M))

#load in UNITED type 2 datasets ----------------------------------------------------------
UNITED_type2 <- create_data(dataset = "united t2d",
                            commonmody = FALSE, id = TRUE)

UNITED <- full_join(UNITED_type1, UNITED_type2)
#Produce characteristics tables ------------------------------------------------------------
## for early-insulin-treated ----------------------------------------------------------------
#create varlist (numeric variables of interest names)
varlist_T1D = c("agedx",
                "agerec", 
                "bmi", 
                "hba1c",
                "durationfinal")
#create varlist_cat (categorical variables of interest names)
varlist_cat_T1D = c("sex", 
                    "pardm", 
                    "insoroha", 
                    "C", 
                    "A", 
                    "T")

#create table for UNITED T1D by MODY status - chosen to have numeric variables displayed as median [IQR]
var_characteristics(varlist = varlist_T1D, 
                    varlist_cat = varlist_cat_T1D, 
                    dataset = UNITED_type1, 
                    numeric_option = "medianIQR", 
                    group = "M",
                    p_value_testing = FALSE,
                    table_name = "Supplementary_Material/Outputs/UNITED_T1D_table")

## for not-early-insulin-treated ----------------------------------------------------------
#create varlist (numeric variables of interest names)
varlist_T2D = c("agedx",
                "agerec", 
                "bmi", 
                "hba1c",
                "durationfinal")
#create varlist_cat (categorical variables of interest names)
varlist_cat_T2D = c("sex", 
                    "pardm", 
                    "insoroha")

#create table for UNITED T1D by MODY status - chosen to have numeric variables displayed as median [IQR]
var_characteristics(varlist = varlist_T2D, 
                    varlist_cat = varlist_cat_T2D, 
                    dataset = UNITED_type2, 
                    numeric_option = "medianIQR", 
                    group = "M",
                    p_value_testing = FALSE,
                    table_name = "Supplementary_Material/Outputs/UNITED_T2D_table")


#Finding the mean MODY prob for M+ AND M- participants
#load predictions
predictions_UNITED_type2 <- readRDS("Model_Predictions/predictions_dataset.UNITED_type2_all_genes_new.rds") %>%
  as.data.frame() %>%
  { rownames(.) <- NULL; . } %>%
  column_to_rownames(var = "id")
predictions_UNITED_type2 <- predictions_UNITED_type2[as.character(c(dataset_UNITED_type2$id)), ]
#bind them
UNITED_type2 <- cbind(UNITED_type2,predictions_UNITED_type2)

UNITED_type2 %>%
  filter(M==0) %>%
  summarise(mean = mean(prob))

UNITED_type2 %>%
  filter(M==1) %>%
  summarise(mean = mean(prob))

UNITED %>%
  summarise(
    median = round(median(durationfinal, na.rm = TRUE),2),
    Q1 = round(quantile(durationfinal, probs = 0.25, na.rm = TRUE),2),
    Q3 = round(quantile(durationfinal, probs = 0.75, na.rm = TRUE),2)
  )