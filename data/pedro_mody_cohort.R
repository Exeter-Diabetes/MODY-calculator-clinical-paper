
# Apply the MODY calculator to everyone in prevalent cohort diagnosed aged 1-35 years

# Investigate time to insulin issues

############################################################################################
  
# Setup
library(tidyverse)
library(aurum)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "test-remote",cprdConf = "~/.aurum.yaml")

analysis = cprd$analysis("dpctn_final")

pedro_mody_cohort <- cohort %>%
  analysis$cached("pedro_mody_cohort", unique_indexes="patid")

pedro_mody_cohort_local <- pedro_mody_cohort %>% collect() %>% mutate(patid=as.character(patid))

############################################################################################

# # Get cohort info
# 
# cohort <- cohort %>% analysis$cached("cohort")
# 
# 
# ############################################################################################
# 
# # Define MODY cohort: patients diagnosed with a current Type 1 or Type 2 diagnosis, diagnosed aged 1-35
# ## Do exclude if missing diagnosis date as invalid, as don't know age of diagnosis
# ## Don't exclude if missing BMI/HbA1c before diagnosis (NB: missing BMI if no BMI in adulthood too)
# 
# pedro_mody_cohort <- cohort %>%
#   filter(dm_diag_age>=1 & dm_diag_age<36 & (diabetes_type=="type 1" | diabetes_type=="type 2" | diabetes_type=="mixed; type 1" | diabetes_type=="mixed; type 2")) %>%
#   mutate(hba1c_mmolmol=ifelse(hba1cdate>=diagnosis_date, hba1c, NA),
#          hba1c=(0.09148*hba1c_mmolmol)+2.152,
#          age_at_bmi=datediff(bmidate, dob)/365.25,
#          bmi=ifelse(bmidate>=diagnosis_date & age_at_bmi>=18, bmi, NA),
#          patient_group=ifelse(dm_diag_age<18, "diag_under_18", "diag_atover_18_ins1"),
#          sex=gender-1,
#          insoroha=ifelse(current_oha==1 | current_insulin==1, 1L, 0L)) %>%
#   select(patid, patient_group, agedx=dm_diag_age, sex, hba1c, bmi, pardm=fh_diabetes, insoroha, agerec=age_at_index) %>%
#   analysis$cached("pedro_mody_cohort", unique_indexes="patid")
# 
# pedro_mody_cohort %>% count()   
# #70258
# 
# pedro_mody_cohort %>% group_by(patient_group) %>% count()
# 
# 
# pedro_mody_cohort_local <- pedro_mody_cohort %>% collect() %>% mutate(patid=as.character(patid))
# 
# save(pedro_mody_cohort_local, file="pedro_mody_cohort.Rda")


############################################################################################

# Looking at mean scores

mody_calc_results <- mody_calc_results %>% analysis$cached("mody_calc_results")

analysis = cprd$analysis("dpctn")
pedro_mody_t1_results <- pedro_mody_t1_results %>% analysis$cached("pedro_mody_t1_results")
pedro_mody_t2_results <- pedro_mody_t2_results %>% analysis$cached("pedro_mody_t2_results")
## Numbers are adjusted MODY probabilities as a proportion out of 1
## Only for those with no missing BMI or HbA1c (n=65172 as per mody_calc_results)
## mean_pardm_0 = assume missing FH=0; mean_pardm_1 = assume missing FH=1

analysis = cprd$analysis("dpctn_final")


#local_data <- mody_calc_results %>% inner_join((pedro_mody_t1_results %>% select(patid, mody_t1_prob=mean_pardm_0)), by="patid") %>% inner_join((pedro_mody_t2_results %>% select(patid, mody_t2_prob=mean_pardm_0)), by="patid") %>% select(patid, which_equation, mody_adj_prob_fh0, mody_t1_prob, mody_t2_prob) %>% mutate(pedro_prob=ifelse(which_equation=="t1", 100*mody_t1_prob, 100*mody_t2_prob)) %>% collect()

local_data_under_30 <- mody_calc_results %>% inner_join((pedro_mody_t1_results %>% select(patid, mody_t1_prob=mean_pardm_0)), by="patid") %>% inner_join((pedro_mody_t2_results %>% select(patid, mody_t2_prob=mean_pardm_0)), by="patid") %>% filter(dm_diag_age<30) %>% select(patid, which_equation, mody_adj_prob_fh0, mody_t1_prob, mody_t2_prob) %>% mutate(pedro_prob=ifelse(which_equation=="t1", 100*mody_t1_prob, 100*mody_t2_prob)) %>% collect()


mean(local_data$mody_adj_prob_fh0)
#8.78%

mean(local_data$pedro_prob)
#3.85%


mean((local_data %>% filter(which_equation=="t1"))$mody_adj_prob_fh0)
#2.40%

mean((local_data %>% filter(which_equation=="t1"))$pedro_prob)
#0.32%


mean((local_data %>% filter(which_equation=="t2"))$mody_adj_prob_fh0)
#13.97%

mean((local_data %>% filter(which_equation=="t2"))$pedro_prob)
#6.71%
