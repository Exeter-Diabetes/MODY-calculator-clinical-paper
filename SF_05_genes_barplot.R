###########################################################################################
#Supplementary Figure 5:
#Stacked barchart of genes in paper datasets

#################################################################################################
#load libraries
library(tidyverse)
#load functions
source("data/create_data.R")

################################################################################################################
#load data 
# MyDiabetes ------------------------------------------------------------------------------------
#load T1D
load("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper/MY_T1D.RData")
#load T2D 
load("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper/MY_T2D.RData")
# LIMIT TO ONLY WHITE ETHNICITY
MY_T1D <- MY_T1D %>%
  filter(ethnicity == "White") %>%
  mutate(M = ifelse(biomark_status == 0, M, NA))

table(MY_T1D$M)
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
dataset.UNITED_p <- create_data(dataset = "united t1d pediatrics", commonmody = FALSE, biomarkers = "reduced")

#need to change M=NA to M=0
dataset.UNITED_p <- dataset.UNITED_p %>%
  filter(!is.na(T))
#checked if worked: should have M=1 (n=10) & M=0 (n=415)
table(dataset.UNITED_p$M)

#T2D
dataset.UNITED_type2p <- dataset.UNITED_p %>%
  filter(tti == "" | tti == "Greater than 12 months")
table(dataset.UNITED_type2p$M)

#T1D
dataset.UNITED_type1p <- dataset.UNITED_p %>%
  filter(!(id %in% dataset.UNITED_type2p$id))
table(dataset.UNITED_type1p$M)

# External dataset -----------------------------------------------------------------------

# Early-insulin-treated
dataset.MYDIABETES_type1 <- MY_T1D %>%
  mutate(study = "MYDIABETES") %>%
  select(MY_ID, study, agerec, agedx, sex, bmi, pardm, insoroha, hba1c, C, A, M, Gene, biomark_status) 

dataset.UNITED_type1p <- dataset.UNITED_type1p %>%
  mutate(study = "UNITED paediatric") %>%
  filter(!is.na(T))
dataset_type1 <- full_join(dataset.MYDIABETES_type1, dataset.UNITED_type1p, by = c("study","agerec", "agedx", "sex", "bmi", "pardm", "insoroha", "hba1c", "C", "A", "M", "Gene")) %>%
  mutate(
    Gene = ifelse(Gene == "", NA, Gene), 
    biomark_status = ifelse(is.na(biomark_status), 
                            ifelse(C == 1 & A == 0, 0, 1),
                            biomark_status)) 
# Not-early-insulin-treated
dataset.MYDIABETES_type2 <- MY_T2D %>%
  mutate(study = "MYDIABETES") %>%
  select(MY_ID, study, agerec, agedx, sex, bmi, pardm, insoroha, hba1c, C, A, M, Gene, biomark_status) 

dataset.UNITED_type2p <- dataset.UNITED_type2p %>%
  mutate(study = "UNITED paediatric") %>%
  select(id, study, agerec, agedx, sex, bmi, pardm, insoroha, hba1c, C, A, M, Gene) 

dataset_type2 <- full_join(dataset.MYDIABETES_type2, dataset.UNITED_type2p, by = c("study","agerec", "agedx", "sex", "bmi", "pardm", "insoroha", "hba1c", "C", "A", "M", "Gene")) %>%
  mutate(
    Gene = ifelse(Gene == "", NA, Gene), 
    biomark_status = ifelse(is.na(biomark_status), 
                            ifelse(C == 1 & A == 0, 0, 1),
                            biomark_status)) 

# case control --------------------------------------------------------------------------------
dataset.cc_type1 <- create_data(dataset = "case-control t1d")
#load in UNITED type 2 dataset
dataset.cc_type2 <- create_data(dataset = "case-control t2d") 

# UNITED --------------------------------------------------------------------------------------
dataset.UNITED_type1 <- create_data(dataset = "united t1d", biomarkers = "full", commonmody = FALSE)

#need to change M=NA to M=0
dataset.UNITED_type1 <- dataset.UNITED_type1 %>%
  mutate(M = ifelse(is.na(M) == TRUE, 0, M))
dataset.UNITED_type2 <- create_data(dataset = "united t2d", commonmody = FALSE)


#--------------------------------------------------------------------------------------
#Make combined dataset
# ---------------------------------------------------------------------------------
  
externalval_type1 <- dataset_type1 %>%
  filter(M==1) %>%
  select(Gene) %>%
  mutate(dataset = "External validation", 
         model = "early insulin")

externalval_type2 <- dataset_type2 %>%
  filter(M == 1) %>%
  select(Gene) %>%
  mutate(dataset = "External validation", 
         model = "not early insulin")

UNITED_type1 <- dataset.UNITED_type1 %>%
  filter(M==1) %>%
  select(Gene) %>%
  mutate(dataset = "Recalibration (UNITED)", 
         model = "early insulin")

UNITED_type2 <- dataset.UNITED_type2 %>%
  filter(M == 1) %>%
  select(Gene) %>%
  mutate(dataset = "Recalibration (UNITED)", 
         model = "not early insulin")

cc_type1 <- dataset.cc_type1 %>%
  filter(typedm != "Type1") %>%
  rename(Gene = typedm) %>%
  select(Gene) %>%
  mutate(dataset = "Training (case-control)", 
         model = "early insulin")

cc_type2 <- dataset.cc_type2 %>%
  filter(typedm != "Type2") %>%
  rename(Gene = typedm) %>%
  select(Gene) %>%
  mutate(dataset = "Training (case-control)", 
         model = "not early insulin")

#Make dataset
genes_data <- rbind(cc_type1, cc_type2, UNITED_type1, UNITED_type2, externalval_type1, externalval_type2)
genes_data$dataset <- factor(genes_data$dataset, levels = c("Training (case-control)", "Recalibration (UNITED)", "External validation"))
genes_data$model <- factor(genes_data$model, levels = c("early insulin", "not early insulin"))




###############################################################################################
#Option 2:
#datasets by genes by model
genes_data <- rbind(cc_type1, cc_type2, UNITED_type1, UNITED_type2, externalval_type1, externalval_type2)
genes_data$dataset <- factor(genes_data$dataset, levels = c("Training (case-control)", "Recalibration (UNITED)", "External validation"))
genes_data$model <- factor(genes_data$model, levels = c("early insulin", "not early insulin"))
genes_data <- genes_data %>%
  group_by(dataset, model) %>%
  count(Gene)


#################################################################################################
#Option 3:
#datasets by genes stacked by model 
#with labels

#this is the one gone for
supfig5 <- ggplot(genes_data, aes(x = Gene, fill = factor(model, labels = c("Early-insulin-treated", "Not-early-insulin-treated")), y = n)) + 
    geom_col(position = "stack") + 
    coord_flip() +
    labs(y = "Count") +
    geom_text(aes(label = n, group = model), position = position_stack(0.5),
              vjust = 0.5) +
    theme_classic() +
    theme(
      strip.background = element_rect(colour="white", fill="white"), 
      strip.text = element_text(size = 12),
      legend.position = "bottom", 
      legend.title = element_blank()
    ) +
    facet_wrap(~dataset, nrow = 3, ncol = 1, scales = "free") 
plot(supfig5)

pdf("figures/supfig5.pdf", width = 6, height = 9)
supfig5
dev.off()
  

#########################################################################################