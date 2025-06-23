#:----------------------------------------------------------------------------
#
#  create_data()
#
#  This file contains a function that engineers the datasets used in this paper
#  into the the right shape 
#
#:----------------------------------------------------------------------------

# ARGUMENTS
# dataset - ("case-control t1d", 
#            "case-control t2d",
#            "united t1d",
#            "united t2d",
#            "united t1d pediatrics") 
#           Which dataset to load and clean
# biomarkers - ("reduced", )
# commonmody - (TRUE (default)/ FALSE) Whether to define MODY genes based on 
#               "common genes" (HNF1A, HNF4A or GCK) or to extend definition 
#               to additional genes
# id - (TRUE / FALSE (default)) Whether to include id numbers in cleaned output


create_data <- function(
    dataset = NULL, 
    biomarkers = "reduced", 
    commonmody = TRUE, 
    id = FALSE) {
  
  
  # FUNCTION USE CHECKS
  
  # checks for 'dataset'
  if(missing(dataset) | is.null(dataset)) {stop("'dataset' needs to be provided.")}
  if(!(tolower(dataset) %in% c("case-control t1d", "case-control t2d", "united t1d", "united t2d", "united t1d pediatrics"))) {stop("'dataset' needs to be: case-control t1d / case-control t2d / united t1d / united t2d / united t1d pediatrics")}
  # checks for 'biomarkers'
  if(!(tolower(biomarkers) %in% c("reduced", "full"))) {stop("'biomarkers' needs to be: reduced / full")}
  # check forr 'commonmody'
  if(!(commonmody %in% c(TRUE, FALSE))) {stop("'commonmody' needs to be: TRUE / FALSE")}
  
  # load libraries
  library(tidyverse)
  
  #### Generate Case-control T1D dataset
  if (tolower(dataset) == "case-control t1d") {
    
    ## load libraries
    library(readr)
    
    ## load dataset
    dataset.case_control <- read_csv("data/mmoct11.csv")
    
    ## select the right patients
    dataset.case_control_type1 <- dataset.case_control %>%
      
      ### make sure outcome is a factor
      mutate(mody = factor(mody)) %>%
      
      ### select those insulin-treated
      filter(t1t2rcgp == 1) %>%
      
      ### select the right variables
      select(mody, sex, bmi, agedx, hba1c, pardm, agerec, typedm, duration) %>%
      
      ### make sure it is a data.frame() for use in R
      as.data.frame()
    
    ## return final dataset
    return(dataset.case_control_type1)
    
  } else if (tolower(dataset) == "case-control t2d") {
    
    ## load libraries
    library(readr)
    
    ## load dataset
    dataset.case_control <- read_csv("data/mmoct11.csv")
    
    ## select the right patients
    dataset.case_control_type2 <- dataset.case_control %>%
      
      ### make sure outcome is a factor
      mutate(mody = factor(mody)) %>%
      
      ### select those non-insulin-treated
      filter(t1t2rcgp == 2) %>%
      
      ### select the right variables
      select(mody, agedx, bmi, hba1c, pardm, agerec, insoroha, sex, typedm, duration) %>%
      
      ### make sure it is a data.frame() for use in R
      as.data.frame()
    
    ## return final dataset
    return(dataset.case_control_type2)
    
  } else if (tolower(dataset) == "united t1d") {
    
    if (tolower(biomarkers) == "reduced") {
      
      ## load libraries
      library(readr)
      
      ## load dataset
      dataset.UNITED <- read_csv("data/UNITEDfull.csv") %>%
        
        ### drop specific Genes that shouldn't be included: keep those without a Gene or those that aren't on this list
        filter(is.na(Gene) | !(Gene %in% c("MIDD", 
                                           "TRMT10A", 
                                           "3243", 
                                           "MDP", 
                                           "MT-ND1", 
                                           "Mitochondrial", 
                                           "WFS1", 
                                           "SLC19A2", 
                                           "ZNF808", 
                                           "INSR", 
                                           "LMNA", 
                                           "PPARG")))
      
      ## select the right patients
      dataset.UNITED_type1 <- dataset.UNITED %>%
        
        ### select those insulin-treated
        filter(t1ort2 == 1) %>%
        
        ### select the correct type of MODY
        filter(commonmody == 1 | 
                 (completepway == 1 & monogenicinctngs == 0 & 
                    is.na(knowncause))) %>%
        
        ### select the correct variables
        select(commonmody, 
               sex, 
               bmi, 
               agedx = agediag, 
               hba1c = hba1cpc, 
               pardm, 
               agerec, 
               UCPCRPosNegFinal, 
               AntibodyFINAL, 
               DNAResult, 
               insoroha, 
               id,
               tti,
               durationfinal) %>%
        
        ### drop patients without patient history for diabetes
        drop_na(pardm) %>%
        
        ### generate the outcome variables
        mutate(mody = DNAResult) %>%
        
        ### drop not needed variables
        select(-commonmody, -DNAResult) %>%
        
        ### drop patients with missing data in the variables we care about
        drop_na(sex, bmi, agedx, hba1c, pardm, agerec) %>%
        
        ### rename variables for biomarkers
        rename(C = UCPCRPosNegFinal, A = AntibodyFINAL, M = mody) %>%
        
        #:--- Lets see it in turns:
        #:---   For M-:
        #:---     - 2 were missing C (C+ as we know A == 0) #:----------------------------------------------------------
      mutate(C = ifelse(M == 0 & A == 0 & is.na(C), 1, C)) %>%
        #:---   For M+:
        #:---     - 1 was C+ and missing A (needs to be A-) #:--------------------------------------------------------
      mutate(A = ifelse(M == 1 & C == 1 & is.na(A), 0, A)) %>%
        #:---     - 1 was missing C and missing A (needs to be C+ and A-)
        mutate(A = ifelse(M == 1 & is.na(C) & is.na(A), 0, A)) %>%
        mutate(C = ifelse(M == 1 & is.na(C) & A == 0, 1, C)) %>%
        #:---   For missing M:
        #:---     - it doesn't matter, there is no one with C+ and A-
        mutate(T = ifelse(C == 0 | A == 1, 1, 0)) %>% # T is 1 if Cn or Ap
        
        ### make sure it is a data.frame() for use in R
        as.data.frame()
      
      
      
      if (!isTRUE(commonmody)) {
        # Find patients with Genes
        # Remove them from the previous dataset to prevent duplication
        # add them back into the dataset with new values
        
        interim_genes <- dataset.UNITED %>%
          
          ### select those insulin-treated
          filter(t1ort2 == 1) %>%
          
          ### keep only patients with Genes
          filter(Gene %in% c("HNF1A", 
                             "HNF4A", 
                             "HNF1B", 
                             "GCK", 
                             "NeuroD1", 
                             "KCNJ11", 
                             "GATA6", 
                             "ABCC8", 
                             "INS", 
                             "RFX6")) %>%
          
          ### select the correct variables
          select(commonmody, 
                 sex,
                 bmi, 
                 agedx = agediag, 
                 hba1c = hba1cpc, 
                 pardm, 
                 agerec, 
                 UCPCRPosNegFinal, 
                 AntibodyFINAL, 
                 DNAResult, 
                 insoroha, 
                 Gene, 
                 id,
                 tti,
                 durationfinal) %>%
          
          ### drop patients without patient history for diabetes
          drop_na(pardm) %>%
          
          ### generate the outcome variables
          mutate(mody = ifelse(!is.na(Gene) 
                               & DNAResult == 1, 1, DNAResult)) %>%
          
          ### drop not needed variables
          select(-commonmody, -DNAResult) %>%
          
          ### drop patients with missing data in the variables we care about
          drop_na(sex, bmi, agedx, hba1c, pardm, agerec) %>%
          
          ### rename variables for biomarkers
          rename(C = UCPCRPosNegFinal, A = AntibodyFINAL, M = mody) %>%
          
          #:--- Lets see it in turns:
          #:---   For M-:
          #:---     - 2 were missing C (C+ as we know A == 0) #:----------------------------------------------------------
        mutate(C = ifelse(M == 0 & A == 0 & is.na(C), 1, C)) %>%
          #:---   For M+:
          #:---     - 1 was C+ and missing A (needs to be A-) #:--------------------------------------------------------
        mutate(A = ifelse(M == 1 & C == 1 & is.na(A), 0, A)) %>%
          #:---     - 1 was missing C and missing A (needs to be C+ and A-)
          mutate(A = ifelse(M == 1 & is.na(C) & is.na(A), 0, A)) %>%
          mutate(C = ifelse(M == 1 & is.na(C) & A == 0, 1, C)) %>%
          #:---   For missing M:
          #:---     - it doesn't matter, there is no one with C+ and A-
          mutate(T = ifelse(C == 0 | A == 1, 1, 0)) %>% # T is 1 if Cn or Ap
          
          ### make sure it is a data.frame() for use in R
          as.data.frame()
        
        # join additional MODY cases
        dataset.UNITED_type1 <- dataset.UNITED_type1 %>%
          
          # remove repeated patients
          filter(!(id %in% interim_genes$id)) %>%
          
          # join them back in
          rbind(
            
            interim_genes %>%
              
              # remove Gene variable
              select(-Gene)
            
          ) %>%
          
          # join Gene variable back in but for everyone
          left_join(
            
            dataset.UNITED %>%
              
              # select id and Gene
              select(id, Gene),
            
            by = c("id")
          )
        
      }
      
      ## return final dataset
      return(
        
        if (isTRUE(id)) {
          
          dataset.UNITED_type1
          
        } else {
          
          dataset.UNITED_type1 %>%
            
            # remove id
            select(-id)
          
        }
        
      )
      
    } else if (tolower(biomarkers) == "full") {
      
      ## load libraries
      library(readr)
      
      ## load dataset
      dataset.UNITED <- read_csv("data/UNITEDfull.csv") %>%
        
        ### drop specific Genes that shouldn't be included: keep those without a Gene or those that aren't on this list
        filter(is.na(Gene) | !(Gene %in% c("MIDD", 
                                           "TRMT10A", 
                                           "3243", 
                                           "MDP", 
                                           "MT-ND1", 
                                           "Mitochondrial", 
                                           "WFS1", 
                                           "SLC19A2", 
                                           "ZNF808", 
                                           "INSR", 
                                           "LMNA", 
                                           "PPARG")))
      
      ## select the right patients
      dataset.UNITED_type1 <- dataset.UNITED %>%
        
        ### select those insulin-treated
        filter(t1ort2 == 1) %>%
        
        ### select the correct type of MODY
        filter(commonmody == 1 | 
                 (completepway == 1 & 
                    monogenicinctngs == 0 & 
                    is.na(knowncause))) %>%
        
        ### select the correct variables
        select(commonmody, 
               sex, 
               bmi, 
               agedx = agediag, 
               hba1c = hba1cpc, 
               pardm, 
               agerec, 
               UCPCRPosNegFinal, 
               AntibodyFINAL, 
               DNAResult, 
               insoroha, 
               ZNT8pos, 
               GADResult, 
               IA2Result, 
               id,
               tti,
               durationfinal) %>%
        
        ### rename Antibody variables
        rename(
          "ZnT8" = "ZNT8pos",
          "GAD" = "GADResult",
          "IA2" = "IA2Result",
        ) %>%
        
        ### change type of column
        mutate(ZnT8 = as.numeric(ZnT8)) %>%
        
        ### drop patients without patient history for diabetes
        drop_na(pardm) %>%
        
        ### generate the outcome variables
        mutate(mody = DNAResult) %>%
        
        ### drop not needed variables
        select(-commonmody, -DNAResult) %>%
        
        ### drop patients with missing data in the variables we care about
        drop_na(sex, bmi, agedx, hba1c, pardm, agerec) %>%
        
        ### edit antibody testing for specific patients
        mutate(
          ZnT8 = ifelse(id %in% c(892, 8002276), 0, ZnT8),
          GAD = ifelse(id %in% c(892, 8002276), 0, GAD),
          IA2 = ifelse(id %in% c(892, 8002276), 0, IA2)
        ) %>%
        
        ### rename variables for biomarkers
        rename(C = UCPCRPosNegFinal, A = AntibodyFINAL, M = mody) %>%
        
        #:--- Lets see it in turns:
        #:---   For M-:
        #:---     - 2 were missing C (C+ as we know A == 0) #:----------------------------------------------------------
      mutate(C = ifelse(M == 0 & A == 0 & is.na(C), 1, C)) %>%
        #:---   For M+:
        #:---     - 1 was C+ and missing A (needs to be A-) #:--------------------------------------------------------
      mutate(A = ifelse(M == 1 & C == 1 & is.na(A), 0, A)) %>%
        #:---     - 1 was missing C and missing A (needs to be C+ and A-)
        mutate(A = ifelse(M == 1 & is.na(C) & is.na(A), 0, A)) %>%
        mutate(C = ifelse(M == 1 & is.na(C) & A == 0, 1, C)) %>%
        #:---   For missing M:
        #:---     - it doesn't matter, there is no one with C+ and A-
        mutate(T = ifelse(C == 0 | A == 1, 1, 0)) %>% # T is 1 if Cn or Ap
        
        ### make sure it is a data.frame() for use in R
        as.data.frame()
      
      if (!isTRUE(commonmody)) {
        # Find patients with Genes
        # Remove them from the previous dataset to prevent duplication
        # add them back into the dataset with new values
        
        interim_genes <- dataset.UNITED %>%
          
          ### select those insulin-treated
          filter(t1ort2 == 1) %>%
          
          ### keep only patients with Genes
          filter(Gene %in% c("HN1A", 
                             "HNF4A", 
                             "HNF1B", 
                             "GCK", 
                             "NeuroD1", 
                             "KCNJ11", 
                             "GATA6", 
                             "ABCC8", 
                             "INS", 
                             "RFX6")) %>%
          
          ### select the correct variables
          select(commonmody, 
                 sex, 
                 bmi, 
                 agedx = agediag, 
                 hba1c = hba1cpc, 
                 pardm, 
                 agerec, 
                 UCPCRPosNegFinal, 
                 AntibodyFINAL, 
                 DNAResult, 
                 insoroha, 
                 ZNT8pos, 
                 GADResult, 
                 IA2Result, 
                 Gene, 
                 id,
                 tti,
                 durationfinal) %>%
          
          ### rename Antibody variables
          rename(
            "ZnT8" = "ZNT8pos",
            "GAD" = "GADResult",
            "IA2" = "IA2Result",
          ) %>%
          
          ### change type of column
          mutate(ZnT8 = as.numeric(ZnT8)) %>%
          
          ### drop patients without patient history for diabetes
          drop_na(pardm) %>%
          
          ### generate the outcome variables
          mutate(mody = ifelse(!is.na(Gene) & DNAResult == 1, 1, DNAResult)) %>%
          
          ### drop not needed variables
          select(-commonmody, -DNAResult) %>%
          
          ### drop patients with missing data in the variables we care about
          drop_na(sex, bmi, agedx, hba1c, pardm, agerec) %>%
          
          ### edit antibody testing for specific patients
          mutate(
            ZnT8 = ifelse(id %in% c(604, 892, 8002276), 0, ZnT8),
            GAD = ifelse(id %in% c(604, 892, 8002276), 0, GAD),
            IA2 = ifelse(id %in% c(604, 892, 8002276), 0, IA2)
          ) %>%
          
          ### rename variables for biomarkers
          rename(C = UCPCRPosNegFinal, A = AntibodyFINAL, M = mody) %>%
          
          #:--- Lets see it in turns:
          #:---   For M-:
          #:---     - 2 were missing C (C+ as we know A == 0) #:----------------------------------------------------------
        mutate(C = ifelse(M == 0 & A == 0 & is.na(C), 1, C)) %>%
          #:---   For M+:
          #:---     - 1 was C+ and missing A (needs to be A-) #:--------------------------------------------------------
        mutate(A = ifelse(M == 1 & C == 1 & is.na(A), 0, A)) %>%
          #:---     - 1 was missing C and missing A (needs to be C+ and A-)
          mutate(A = ifelse(M == 1 & is.na(C) & is.na(A), 0, A)) %>%
          mutate(C = ifelse(M == 1 & is.na(C) & A == 0, 1, C)) %>%
          #:---   For missing M:
          #:---     - it doesn't matter, there is no one with C+ and A-
          mutate(T = ifelse(C == 0 | A == 1, 1, 0)) %>% # T is 1 if Cn or Ap
          
          ### make sure it is a data.frame() for use in R
          as.data.frame()
        
        
        
        # join additional MODY cases
        dataset.UNITED_type1 <- dataset.UNITED_type1 %>%
          
          # remove repeated patients
          filter(!(id %in% interim_genes$id)) %>%
          
          # join them back in
          rbind(
            
            interim_genes %>%
              
              # remove Gene variable
              select(-Gene)
            
          ) %>%
          
          # join Gene variable back in but for everyone
          left_join(
            
            dataset.UNITED %>%
              
              # select id and Gene
              select(id, Gene),
            
            by = c("id")
          )
        
      }
      
      ## return final dataset
      return(
        
        
        if (isTRUE(id)) {
          
          dataset.UNITED_type1
          
        } else {
          
          dataset.UNITED_type1 %>%
            
            # remove id
            select(-id)
          
        }
        
      )
      
    }
    
    
  } else if (tolower(dataset) == "united t2d") {
    
    ## load libraries
    library(readr)
    
    ## load dataset
    dataset.UNITED <- read_csv("data/UNITEDfull.csv") %>%
      
      ### drop specific Genes that shouldn't be included: keep those without a Gene or those that aren't on this list
      filter(is.na(Gene) | !(Gene %in% c("MIDD", 
                                         "TRMT10A", 
                                         "3243", 
                                         "MDP", 
                                         "MT-ND1", 
                                         "Mitochondrial", 
                                         "WFS1", 
                                         "SLC19A2", 
                                         "ZNF808", 
                                         "INSR", 
                                         "LMNA", 
                                         "PPARG")))
    
    
    ## select the right patients
    dataset.UNITED_type2 <- dataset.UNITED %>%
      
      ### select those non-insulin-treated
      filter(t1ort2 == 2) %>%
      
      ### select the correct type of MODY
      filter(commonmody == 1 | (completepway == 1 & monogenicinctngs == 0 & is.na(knowncause))) %>%
      
      ### select the correct variables
      select(commonmody, 
             sex, 
             bmi, 
             agedx = agediag, 
             insoroha, 
             hba1c = hba1cpc, 
             pardm, 
             agerec, 
             UCPCRPosNegFinal, 
             AntibodyFINAL, 
             DNAResult, 
             id,
             tti,
             durationfinal) %>%
      
      ### generate the outcome variables
      mutate(mody = DNAResult) %>%
      
      ### remove variables not needed
      select(-commonmody, -DNAResult) %>%
      
      ### rename variables
      rename(M = mody) %>%
      
      ### select the variables we needed
      select(c("sex", 
               "bmi", 
               "agedx", 
               "insoroha", 
               "hba1c", 
               "pardm", 
               "agerec", 
               "M", 
               "id", 
               "tti",
               "durationfinal",
               "UCPCRPosNegFinal", 
               "AntibodyFINAL" )) %>%
      
      ### drop patients with missing data in the variables we care about
      drop_na(c("sex", 
                "bmi", 
                "agedx", 
                "insoroha", 
                "hba1c", 
                "pardm", 
                "agerec", 
                "M")) %>%
      
      ### make sure it is a data.frame() for use in R
      as.data.frame() %>%
      
      ### change specific patient
      mutate(M = ifelse(id %in% c(8002013), 1, M))
    
    if (!isTRUE(commonmody)) {
      # Find patients with Genes
      # Remove them from the previous dataset to prevent duplication
      # add them back into the dataset with new values
      
      interim_genes <- dataset.UNITED %>%
        
        ### select those insulin-treated
        filter(t1ort2 == 2) %>%
        
        ### keep only patients with Genes
        filter(Gene %in% c("HNF1A", 
                           "HNF4A", 
                           "HNF1B", 
                           "GCK", 
                           "NeuroD1", 
                           "KCNJ11", 
                           "GATA6", 
                           "ABCC8", 
                           "INS", 
                           "RFX6")) %>%
        
        ### select the correct variables
        select(commonmody, 
               sex, 
               bmi, 
               agedx = agediag, 
               insoroha, 
               hba1c = hba1cpc, 
               pardm, 
               agerec, 
               UCPCRPosNegFinal, 
               AntibodyFINAL, 
               DNAResult, 
               Gene, 
               id, 
               tti,
               durationfinal) %>%
        
        ### change specific patient
        mutate(DNAResult = ifelse(id %in% c(8002013), 1, DNAResult)) %>%
        ### generate the outcome variables
        mutate(mody = ifelse(!is.na(Gene) & DNAResult == 1, 1, DNAResult)) %>%
        
        ### remove variables not needed
        select(-commonmody, -DNAResult) %>%
        
        ### rename variables
        rename(M = mody) %>%
        
        ### select the variables we needed
        select(c("sex", 
                 "bmi", 
                 "agedx", 
                 "insoroha", 
                 "hba1c", 
                 "pardm", 
                 "agerec", 
                 "M", 
                 "Gene", 
                 "id", 
                 "tti",
                 "durationfinal",
                 "UCPCRPosNegFinal", 
                 "AntibodyFINAL")) %>%
        
        ### drop patients with missing data in the variables we care about
        drop_na(c("sex", 
                  "bmi", 
                  "agedx", 
                  "insoroha", 
                  "hba1c", 
                  "pardm", 
                  "agerec", 
                  "M")) %>%
        
        ### make sure it is a data.frame() for use in R
        as.data.frame()
      
      
      # join additional MODY cases
      dataset.UNITED_type2 <- dataset.UNITED_type2 %>%
        
        # remove repeated patients
        filter(!(id %in% interim_genes$id)) %>%
        
        # join them back in
        rbind(
          
          interim_genes %>%
            
            # remove Gene variable
            select(-Gene)
          
        ) %>%
        
        # join Gene variable back in but for everyone
        left_join(
          
          dataset.UNITED %>%
            
            # select id and Gene
            select(id, Gene),
          
          by = c("id")
        )
      
    }
    
    ## return final dataset
    return(
      
      
      if (isTRUE(id)) {
        
        dataset.UNITED_type2 %>%
          
        filter(!(id %in% c(540)))
        
      } else {
        
        dataset.UNITED_type2 %>%
          
          filter(!(id %in% c(540))) %>%
          
          # remove id
          select(-id)
        
      }
      
    )
    
  } else if (tolower(dataset) == "united t1d pediatrics") {
    
    if (tolower(biomarkers) == "reduced") {
      
      ## load libraries
      library(readr)
      library(haven)
      
      ## load dataset
      dataset.UNITEDfull <- read_dta("data/UNITEDfull.dta") %>%
        
        ### drop specific Genes that shouldn't be included: keep those without a Gene or those that aren't on this list
        filter(is.na(Gene) | !(Gene %in% c("MIDD", 
                                           "TRMT10A", 
                                           "3243", 
                                           "MDP", 
                                           "MT-ND1", 
                                           "Mitochondrial", 
                                           "WFS1", 
                                           "SLC19A2", 
                                           "ZNF808", 
                                           "INSR", 
                                           "LMNA", 
                                           "PPARG")))
      
      dataset.UNITED <- read_csv("data/UNITEDfull.csv") %>%
        
        ### drop specific Genes that shouldn't be included: keep those without a Gene or those that aren't on this list
        filter(is.na(Gene) | !(Gene %in% c("MIDD", 
                                           "TRMT10A", 
                                           "3243", 
                                           "MDP", 
                                           "MT-ND1", 
                                           "Mitochondrial", 
                                           "WFS1", 
                                           "SLC19A2", 
                                           "ZNF808", 
                                           "INSR", 
                                           "LMNA", 
                                           "PPARG")))
      
      
      ## select the right patients
      dataset.UNITED_young_type1 <- dataset.UNITEDfull %>% 
        
        ### remove patients used in the other analysis
        filter(!(id %in% dataset.UNITED$id)) %>% 
        filter(swpaed != 0) %>%

        ### select the correct variables
        select(sex, 
               bmi, 
               agedx = agediag, 
               hba1c = hba1cpc, 
               pardm, 
               agerec, 
               UCPCRPosNegFinal, 
               AntibodyFINAL, 
               DNAResult, 
               id, 
               Gene,
               tti,
               insoroha,
               durationfinal) %>%
        
        ### drop patients without patient history for diabetes
        drop_na(pardm) %>%
        
        ### generate the outcome variables
        mutate(mody = ifelse(is.na(DNAResult), 
                             NA, 
                             ifelse(tolower(Gene) %in% c("hnf1a", 
                                                         "hnf4a", 
                                                         "gck") & 
                                      DNAResult == 1, 1, 0))) %>%
        
        ### drop not needed variables
        select(-DNAResult, -Gene) %>%
        
        ### drop patients with missing data in the variables we care about
        drop_na(sex, 
                bmi, 
                agedx, 
                hba1c, 
                pardm, 
                agerec) %>%
        
        ### rename variables for biomarkers
        rename(C = UCPCRPosNegFinal, A = AntibodyFINAL, M = mody) %>%
        
        #:--- Lets see it in turns:
        #:---   For M-:
        #:---     - 2 were missing C (C+ as we know A == 0) #:----------------------------------------------------------
      mutate(C = ifelse(M == 0 & A == 0 & is.na(C), 1, C)) %>%
        #:---   For M+:
        #:---     - 1 was C+ and missing A (needs to be A-) #:--------------------------------------------------------
      mutate(A = ifelse(M == 1 & C == 1 & is.na(A), 0, A)) %>%
        #:---     - 1 was missing C and missing A (needs to be C+ and A-)
        mutate(A = ifelse(M == 1 & is.na(C) & is.na(A), 0, A)) %>%
        mutate(C = ifelse(M == 1 & is.na(C) & A == 0, 1, C)) %>%
        mutate(A = ifelse(M == 0 & C == 1 & is.na(A), 0, A)) %>%
        mutate(C = ifelse(!is.na(M), 1, C)) %>%
        mutate(A = ifelse(!is.na(M), 0, A)) %>%
        mutate(A = ifelse(id %in% c(946), 1, A)) %>%      # patient with wrong antibody status coded (IA2 was positive)
        #:---   For missing M:
        #:---     - it doesn't matter, there is no one with C+ and A-
        mutate(T = ifelse(C == 0 | A == 1, 1, 0)) %>% # T is 1 if Cn or Ap
        ### make sure it is a data.frame() for use in R
        as.data.frame()
      
      
      
      if (!isTRUE(commonmody)) {
        # Find patients with Genes
        # Remove them from the previous dataset to prevent duplication
        # add them back into the dataset with new values
        
        interim_genes <- dataset.UNITEDfull %>%
          
          ### remove patients used in the other analysis
          filter(!(id %in% dataset.UNITED$id)) %>% 
          filter(swpaed != 0) %>%
          
          ### keep only patients with Genes
          filter(Gene %in% c("HN1A", 
                             "HNF4A", 
                             "HNF1B", 
                             "GCK", 
                             "NeuroD1", 
                             "KCNJ11", 
                             "GATA6", 
                             "ABCC8", 
                             "INS", 
                             "RFX6")) %>%
          
          ### select the correct variables
          select(sex, 
                 bmi, 
                 agedx = agediag, 
                 hba1c = hba1cpc, 
                 pardm, 
                 agerec, 
                 UCPCRPosNegFinal, 
                 AntibodyFINAL, 
                 DNAResult, 
                 Gene, 
                 id,
                 tti,
                 insoroha,
                 durationfinal) %>%
          
          ### drop patients without patient history for diabetes
          drop_na(pardm) %>%
          
          ### generate the outcome variables
          mutate(mody = ifelse(!is.na(Gene) & 
                                 DNAResult == 1, 1, DNAResult)) %>%
          
          ### drop not needed variables
          select(-DNAResult) %>%
          
          ### drop patients with missing data in the variables we care about
          drop_na(sex, 
                  bmi, 
                  agedx, 
                  hba1c, 
                  pardm, 
                  agerec) %>%
          
          ### rename variables for biomarkers
          rename(C = UCPCRPosNegFinal, A = AntibodyFINAL, M = mody) %>%
          
          #:--- Lets see it in turns:
          #:---   For M-:
          #:---     - 2 were missing C (C+ as we know A == 0) #:----------------------------------------------------------
        mutate(C = ifelse(M == 0 & A == 0 & is.na(C), 1, C)) %>%
          #:---   For M+:
          #:---     - 1 was C+ and missing A (needs to be A-) #:--------------------------------------------------------
        mutate(A = ifelse(M == 1 & C == 1 & is.na(A), 0, A)) %>%
          #:---     - 1 was missing C and missing A (needs to be C+ and A-)
          mutate(A = ifelse(M == 1 & is.na(C) & is.na(A), 0, A)) %>%
          mutate(C = ifelse(M == 1 & is.na(C) & A == 0, 1, C)) %>%
          mutate(A = ifelse(M == 0 & C == 1 & is.na(A), 0, A)) %>%
          #:---   For missing M:
          #:---     - it doesn't matter, there is no one with C+ and A-
          mutate(T = ifelse(C == 0 | A == 1, 1, 0)) %>% # T is 1 if Cn or Ap
          mutate(A = ifelse(id %in% c(946), 1, A)) %>%      # patient with wrong antibody status coded (IA2 was positive)
          
          ### make sure it is a data.frame() for use in R
          as.data.frame()
        
        # join additional MODY cases
        dataset.UNITED_young_type1 <- dataset.UNITED_young_type1 %>%
          
          # remove repeated patients
          filter(!(id %in% interim_genes$id)) %>%
          
          # join them back in
          rbind(
            
            interim_genes %>%
              
              # remove Gene variable
              select(-Gene)
            
          ) %>%
          
          # join Gene variable back in but for everyone
          left_join(
            
            dataset.UNITEDfull %>%
              
              # select id and Gene
              select(id, Gene),
            
            by = c("id")
          )
        
      }
      
      ## return final dataset
      return(
        
        
        if (isTRUE(id)) {
          
          dataset.UNITED_young_type1 %>%
            
            # drop those with T = NA
            drop_na(T)
          
        } else {
          
          dataset.UNITED_young_type1 %>%
            
            # drop those with T = NA
            drop_na(T) %>%
            
            # remove id
            select(-id)
          
        }
      )
      
    } else if (tolower(biomarkers) == "full") {
      
      ## load libraries
      library(readr)
      library(haven)
      
      ## load dataset
      dataset.UNITEDfull <- read_dta("data/UNITEDfull.dta") %>%
        
        ### drop specific Genes that shouldn't be included: keep those without a Gene or those that aren't on this list
        filter(is.na(Gene) | !(Gene %in% c("MIDD", 
                                           "TRMT10A", 
                                           "3243", 
                                           "MDP", 
                                           "MT-ND1", 
                                           "Mitochondrial", 
                                           "WFS1", 
                                           "SLC19A2", 
                                           "ZNF808", 
                                           "INSR", 
                                           "LMNA", 
                                           "PPARG")))
      
      dataset.UNITED <- read_csv("data/UNITEDfull.csv") %>%
        
        ### drop specific Genes that shouldn't be included: keep those without a Gene or those that aren't on this list
        filter(is.na(Gene) | !(Gene %in% c("MIDD", 
                                           "TRMT10A", 
                                           "3243", 
                                           "MDP", 
                                           "MT-ND1", 
                                           "Mitochondrial", 
                                           "WFS1", 
                                           "SLC19A2", 
                                           "ZNF808", 
                                           "INSR", 
                                           "LMNA", 
                                           "PPARG")))
      
      
      ## select the right patients
      dataset.UNITED_young_type1 <- dataset.UNITEDfull %>% 
        
        ### remove patients used in the other analysis
        filter(!(id %in% dataset.UNITED$id)) %>% 
        filter(swpaed != 0) %>%
        
        ### select the correct variables
        select(sex, 
               bmi, 
               agedx = agediag, 
               hba1c = hba1cpc, 
               pardm, 
               agerec, 
               UCPCRPosNegFinal, 
               AntibodyFINAL, 
               DNAResult, 
               id, 
               Gene,
               ZNT8pos, 
               GADResult, 
               IA2Result, 
               tti,
               durationfinal) %>%
        
        ### rename Antibody variables
        rename(
          "ZnT8" = "ZNT8pos",
          "GAD" = "GADResult",
          "IA2" = "IA2Result",
        ) %>%
        
        ### change type of column
        mutate(ZnT8 = as.numeric(ZnT8)) %>%
        
        ### drop patients without patient history for diabetes
        drop_na(pardm) %>%
        
        ### generate the outcome variables
        mutate(mody = ifelse(is.na(DNAResult), 
                             NA, 
                             ifelse(tolower(Gene) %in% c("hnf1a", 
                                                         "hnf4a", 
                                                         "gck"), 
                                    1, 
                                    0))) %>%
        
        ### drop not needed variables
        select(-DNAResult, -Gene) %>%
        
        ### drop patients with missing data in the variables we care about
        drop_na(sex, 
                bmi, 
                agedx, 
                hba1c, 
                pardm, 
                agerec) %>%
        
        ### rename variables for biomarkers
        rename(C = UCPCRPosNegFinal, A = AntibodyFINAL, M = mody) %>%
        
        #:--- Lets see it in turns:
        #:---   For M-:
        #:---     - 2 were missing C (C+ as we know A == 0) #:----------------------------------------------------------
      mutate(C = ifelse(M == 0 & A == 0 & is.na(C), 1, C)) %>%
        #:---   For M+:
        #:---     - 1 was C+ and missing A (needs to be A-) #:--------------------------------------------------------
      mutate(A = ifelse(M == 1 & C == 1 & is.na(A), 0, A)) %>%
        #:---     - 1 was missing C and missing A (needs to be C+ and A-)
        mutate(A = ifelse(M == 1 & is.na(C) & is.na(A), 0, A)) %>%
        mutate(C = ifelse(M == 1 & is.na(C) & A == 0, 1, C)) %>%
        mutate(A = ifelse(M == 0 & C == 1 & is.na(A), 0, A)) %>%
        mutate(C = ifelse(!is.na(M), 1, C)) %>%
        mutate(A = ifelse(!is.na(M), 0, A)) %>%
        #:---   For missing M:
        #:---     - it doesn't matter, there is no one with C+ and A-
        mutate(T = ifelse(C == 0 | A == 1, 1, 0)) %>% # T is 1 if Cn or Ap
        mutate(A = ifelse(id %in% c(946), 1, A)) %>%      # patient with wrong antibody status coded (IA2 was positive)
        
        ### make sure it is a data.frame() for use in R
        as.data.frame()
      
      
      
      if (!isTRUE(commonmody)) {
        # Find patients with Genes
        # Remove them from the previous dataset to prevent duplication
        # add them back into the dataset with new values
        
        interim_genes <- dataset.UNITEDfull %>%
          
          ### remove patients used in the other analysis
          filter(!(id %in% dataset.UNITED$id)) %>% 
          filter(swpaed != 0) %>%
          
          ### keep only patients with Genes
          filter(Gene %in% c("HN1A", 
                             "HNF4A", 
                             "HNF1B", 
                             "GCK", 
                             "NeuroD1", 
                             "KCNJ11", 
                             "GATA6", 
                             "ABCC8", 
                             "INS", 
                             "RFX6")) %>%
          
          ### select the correct variables
          select(sex, 
                 bmi, 
                 agedx = agediag, 
                 hba1c = hba1cpc, 
                 pardm, 
                 agerec, 
                 UCPCRPosNegFinal, 
                 AntibodyFINAL, 
                 DNAResult, 
                 Gene, 
                 id,
                 ZNT8pos, 
                 GADResult, 
                 IA2Result, 
                 tti,
                 durationfinal) %>%
          
          ### rename Antibody variables
          rename(
            "ZnT8" = "ZNT8pos",
            "GAD" = "GADResult",
            "IA2" = "IA2Result",
          ) %>%
          
          ### change type of column
          mutate(ZnT8 = as.numeric(ZnT8)) %>%
          
          ### drop patients without patient history for diabetes
          drop_na(pardm) %>%
          
          ### generate the outcome variables
          mutate(mody = ifelse(!is.na(Gene), 1, DNAResult)) %>%
          
          ### drop not needed variables
          select(-DNAResult) %>%
          
          ### drop patients with missing data in the variables we care about
          drop_na(sex, 
                  bmi, 
                  agedx, 
                  hba1c, 
                  pardm, 
                  agerec) %>%
          
          ### rename variables for biomarkers
          rename(C = UCPCRPosNegFinal, A = AntibodyFINAL, M = mody) %>%
          
          #:--- Lets see it in turns:
          #:---   For M-:
          #:---     - 2 were missing C (C+ as we know A == 0) #:----------------------------------------------------------
        mutate(C = ifelse(M == 0 & A == 0 & is.na(C), 1, C)) %>%
          #:---   For M+:
          #:---     - 1 was C+ and missing A (needs to be A-) #:--------------------------------------------------------
        mutate(A = ifelse(M == 1 & C == 1 & is.na(A), 0, A)) %>%
          #:---     - 1 was missing C and missing A (needs to be C+ and A-)
          mutate(A = ifelse(M == 1 & is.na(C) & is.na(A), 0, A)) %>%
          mutate(C = ifelse(M == 1 & is.na(C) & A == 0, 1, C)) %>%
          mutate(A = ifelse(M == 0 & C == 1 & is.na(A), 0, A)) %>%
          #:---   For missing M:
          #:---     - it doesn't matter, there is no one with C+ and A-
          mutate(T = ifelse(C == 0 | A == 1, 1, 0)) %>% # T is 1 if Cn or Ap
          mutate(A = ifelse(id %in% c(946), 1, A)) %>%      # patient with wrong antibody status coded (IA2 was positive)
          
          ### make sure it is a data.frame() for use in R
          as.data.frame()
        
        # join additional MODY cases
        dataset.UNITED_young_type1 <- dataset.UNITED_young_type1 %>%
          
          # remove repeated patients
          filter(!(id %in% interim_genes$id)) %>%
          
          # join them back in
          rbind(
            
            interim_genes %>%
              
              # remove Gene variable
              select(-Gene)
            
          ) %>%
          
          # join Gene variable back in but for everyone
          left_join(
            
            dataset.UNITEDfull %>%
              
              # select id and Gene
              select(id, Gene),
            
            by = c("id")
          )
        
      }
      
      ## return final dataset
      return(
        
        
        
        if (isTRUE(id)) {
          
          dataset.UNITED_young_type1 %>%
            
            # drop those with T = NA
            drop_na(T)
          
        } else {
          
          dataset.UNITED_young_type1 %>%
            
            # drop those with T = NA
            drop_na(T) %>%
            
            # remove id
            select(-id)
          
        }
        
      )
      
      
      
    }
    
    
  }
  
  
  
  # referral t1d
  
  ## need to load data from the github with Julieanne
  ## need to load script
  ## format the rest
  
  # referral t2d
  
  ## need to load data from the github with Julieanne
  ## need to load script
  ## format the rest
  
  
  
  
  
}
