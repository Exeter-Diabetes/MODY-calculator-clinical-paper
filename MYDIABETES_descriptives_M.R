##MODY in MYDiabetes -------------------------------------------------------------------------------------------------------------------------

##read in libraries --------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(readxl)
library(writexl)

##data preparation ----------------------------------------------------------------------------------------------------------------------
### read in data ------------------------------------------------------------------------------------------------------------------------
MYDIABETES <- read_xlsx("MYDIABETES.xlsx")


### overall summary of data --------------------------------------------------------------------------------------------------------------
summary(MYDIABETES)
### explanation of variable names -------------------------------------------------------------------------------------------------------
#Age = Age at recruitment
#eth = ethnicity
#anonid = ID
#smactualresult2 = MODY status
#smGene = MODY gene


### rename variables for understandability and to fit modelling method ------------------------------------------------------------------
MYDIABETES <- MYDIABETES %>%
  rename(
    MODY = smactualresult2,
    ethnicity = eth, 
    Gene = smGene,
    MY_ID = anonid, 
    agedx = age_diagnosis,
    agerec = Age, 
    Cpeptide = `cpeptide (pmol/L)`, 
    hba1c = `hba1c (mmol/mol)`
  )


### Make recoded version of Time_ins_from_diag ------------------------------------------------------------------------------------------------------
MYDIABETES <- MYDIABETES %>%
  mutate(insfromdiag = ifelse(!is.na(Time_Ins_from_diag), 
                              Time_Ins_from_diag,
                              ifelse(is.na(Initial_Treatment) & is.na(Current_Treatment),
                                     NA,
                                     ifelse(((!is.na(Initial_Treatment) & Initial_Treatment == "Tablets") | (!is.na(Initial_Treatment) & Initial_Treatment == "Diet only")) & ((!is.na(Current_Treatment) & Current_Treatment == "Tablets") | (!is.na(Current_Treatment) & Current_Treatment == "Diet Only")),
                                            "Never on insulin",
                                            ifelse(((!is.na(Initial_Treatment) & Initial_Treatment == "Insulin") | (!is.na(Initial_Treatment) & Initial_Treatment == "Tablets & Insulin")) & ((!is.na(Current_Treatment) & Current_Treatment == "Insulin") | (!is.na(Current_Treatment) &Current_Treatment == "Tablets & Insulin")),
                                                   "Initial and Current treat is ins unknown time",
                                                   ifelse((!is.na(Initial_Treatment) & Initial_Treatment == "Insulin") | (!is.na(Initial_Treatment) & Initial_Treatment == "Tablets & Insulin"),
                                                          "Initial treat only is ins unknown time",
                                                          ifelse((!is.na(Current_Treatment) & Current_Treatment == "Insulin") | (!is.na(Current_Treatment) & Current_Treatment == "Tablets & Insulin"),
                                                                 "Current treat only is ins unknown time",
                                                                 "Not currently on insulin unknown time or initial"
                                                                 )
                                                          )
                                                   )
                                            )
                                     )
                              )
         )
         
table(MYDIABETES$insfromdiag, useNA = "ifany")



### make new variables -------------------------------------------------------------------------------------------------------------------
MYDIABETES <- MYDIABETES %>%
  mutate(
    #pardm = at least 1 parent affected with diabetes
    pardm = ifelse(FH_Mother == "Yes" | FH_Father == "Yes", 
                   "Yes",
                   "No"), 
    #this decides which model arm to go down for calcualting MODY probability
    #if go onto insulin within 6 (12) months - early-insulin-treated model (or "T1D")
    #if don't go onto insulin or after 6 (12) months - not-early-insulin-treated model (or "T2D")
    mody_model = ifelse(insfromdiag == "<12 months" | insfromdiag == "Immediately" | insfromdiag == "Current treat only is ins unknown time" | insfromdiag == "Initial and Current treat is ins unknown time", 
                        "T1D", 
                        "T2D"), 
    #defines C (cpeptide status) as 1 if >= 200 pmol/l and 0 if less than 200 pmol/l
    C = ifelse(Cpeptide >= 80, 1, 0), 
    #defines A (antibody status) as 1 if 1+ islet autoantibody positive and 0 if all negative
    A = ifelse(is.na(GAD_PN) & is.na(IA2_PN) & is.na(ZnT8_PN0),
               NA,
               ifelse((!is.na(ZnT8_PN0) & ZnT8_PN0 == "Positive") | (!is.na(GAD_PN) & GAD_PN == "Positive") | (!is.na(IA2_PN) & IA2_PN == "Positive"),
                      1,
                      0)),
    #this variable denotes which biomarker status (based on cpeptide (C) and islet autibodies(A)) individual has
    #0 is if Cpep >= 200pmol/l AND A == 0 (all negative), 1 is if Cpep < 200 pmol/l OR A ==1 (1+ positive antibody)
    biomark_status = ifelse(C == 1 & A == 0, 0, 1), 
    #if diagnosed under age of 30
    agediag30 = ifelse(agedx < 30, "Yes", "No"),
    #if diagnosed under age of 35
    agediag35 = ifelse(agedx < 35, "Yes", "No"),
    #convert hba1c mmol/mol into %
    hba1c_perc = (0.0915*hba1c)+2.15,
    #clean up M
    M = ifelse(MODY == "Mutation" | MODY == "VUS", 1, 0),
    #Make MODY variable based on 3 genes only
    M3 = ifelse(is.na(Gene),
                0,
                ifelse(Gene == "GCK" | Gene == "HNF1A" | Gene == "HNF4A", 
                       1, 
                       0)),
    #if went to insulin within 12 months
    insin12 = ifelse((!is.na(Time_Ins_from_diag) & Time_Ins_from_diag == "<12 months") | (!is.na(Time_Ins_from_diag) & Time_Ins_from_diag == "Immediately"), 
                     "Yes", 
                     "No"), 
    insoroha = ifelse(Current_Treatment == "Insulin" | Current_Treatment == "Tablets" | Current_Treatment == "Tablets & Insulin", "Yes", "No" )
  )


MYDIABETES %>%
  filter(is.na(Time_Ins_from_diag) & (Initial_Treatment == "Insulin" | Initial_Treatment == "Tablets & Insulin" | Current_Treatment == "Insulin" | Current_Treatment == "Tablets & Insulin")) %>%
  select(MY_ID, ethnicity, Gene, agedx, agerec, Cpeptide, C, Time_Ins_from_diag, Initial_Treatment, Current_Treatment, A, GAD_PN, GAD, IA2_PN, IA2, ZnT8_PN0, ZnT8) %>%
  print(n=76)
MYDIABETES %>%
  filter(is.na(insfromdiag)) %>%
  select(MY_ID, ethnicity, Gene, agedx, agerec, Cpeptide, C, Time_Ins_from_diag, Initial_Treatment, Current_Treatment, A, GAD_PN, GAD, IA2_PN, IA2, ZnT8_PN0, ZnT8) %>%
  print(n=76)

MYDIABETES %>%
  filter(insfromdiag == "Not currently on insulin unknown time or initial") %>%
  select(MY_ID, ethnicity, Gene, agedx, agerec, Cpeptide, C, Time_Ins_from_diag, Initial_Treatment, Current_Treatment, A, GAD_PN, GAD, IA2_PN, IA2, ZnT8_PN0, ZnT8) %>%
  print(n=76)

MYDIABETES %>%
  filter(is.na(Time_Ins_from_diag) & (Initial_Treatment == "Insulin" | Initial_Treatment == "Tablets & Insulin" | Current_Treatment == "Insulin" | Current_Treatment == "Tablets & Insulin")) %>%
  select(MY_ID, ethnicity, Gene, agedx, agerec, Cpeptide, C, Time_Ins_from_diag, Initial_Treatment, Current_Treatment, A, GAD_PN, GAD, IA2_PN, IA2, ZnT8_PN0, ZnT8) %>%
  print(n=76)

agedxmissing <- MYDIABETES %>%
  filter(is.na(agedx)) %>%
  select(MY_ID, ethnicity, Gene, agedx, agerec, Cpeptide, C, Time_Ins_from_diag, Initial_Treatment, Current_Treatment, A, GAD_PN, GAD, IA2_PN, IA2, ZnT8_PN0, ZnT8) %>%
  print(n=76)

currentmissing <- MYDIABETES %>%
  filter(is.na(Current_Treatment)) %>%
  select(MY_ID, ethnicity, Gene, agedx, agerec, Cpeptide, C, Time_Ins_from_diag, Initial_Treatment, Current_Treatment, A, GAD_PN, GAD, IA2_PN, IA2, ZnT8_PN0, ZnT8) %>%
  print(n=76)

write_xlsx(agedxmissing,"agedxmissing_table.xlsx")
write_xlsx(currentmissing,"currentmissing_table.xlsx")


save(MYDIABETES, file = "MYDIABETES.RData")

#Load mydiabetes data
load("MYDIABETES.RData")
## whole data characteristics by MODY status ---------------------------------------------------------------------------------------------
source("var_characteristics.R")

#create varlist (numeric variables of interest names)
