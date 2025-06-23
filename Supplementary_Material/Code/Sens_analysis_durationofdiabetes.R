########################################################################################
#Duration of diabetes in datasets
#######################################################################################
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
CC <- full_join(cc_type1, cc_type2)
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
                            commonmody = FALSE)

UNITED <- full_join(UNITED_type1, UNITED_type2)
# MyDiabetes ------------------------------------------------------------------------------------
#load T1D
load("Data/MY_T1D.RData")
#load T2D 
load("Data/MY_T2D.RData")
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
UNITED_p <- create_data(dataset = "united t1d pediatrics", 
                        commonmody = FALSE, 
                        biomarkers = "reduced",
                        id = TRUE)

#need to change M=NA to M=0
UNITED_p <- UNITED_p %>%
  filter(!is.na(T))
#checked if worked: should have M=1 (n=10) & M=0 (n=415)
table(UNITED_p$M)

#T2D
UNITED_type2p <- UNITED_p %>%
  filter(tti == "" | tti == "Greater than 12 months")
table(UNITED_type2p$M)

#T1D
UNITED_type1p <- UNITED_p %>%
  filter(!(id %in% UNITED_type2p$id))
table(UNITED_type1p$M)

# External dataset -----------------------------------------------------------------------
## Early-insulin-treated
MYDIABETES_type1 <- MY_T1D %>%
  mutate(study = "MYDIABETES") %>%
  select(MY_ID, study, agerec, agedx, sex, bmi, pardm, insoroha, hba1c, C, A, M, Gene, biomark_status) 

UNITED_type1p <- UNITED_type1p %>%
  mutate(study = "UNITED paediatric") %>%
  filter(!is.na(T))
dataset_type1 <- full_join(MYDIABETES_type1, 
                           UNITED_type1p, 
                           by = c("study",
                                  "agerec", 
                                  "agedx", 
                                  "sex", 
                                  "bmi", 
                                  "pardm", 
                                  "insoroha", 
                                  "hba1c", 
                                  "C", 
                                  "A", 
                                  "M", 
                                  "Gene")) %>%
  mutate(
    Gene = ifelse(Gene == "", NA, Gene), 
    biomark_status = ifelse(is.na(biomark_status), 
                            ifelse(C == 1 & A == 0, 0, 1),
                            biomark_status)) 
# Not-early-insulin-treated
MYDIABETES_type2 <- MY_T2D %>%
  mutate(study = "MYDIABETES") %>%
  select(MY_ID, study, agerec, agedx, sex, bmi, pardm, insoroha, hba1c, C, A, M, Gene, biomark_status) 

UNITED_type2p <- UNITED_type2p %>%
  mutate(study = "UNITED paediatric") %>%
  select(id, study, agerec, agedx, sex, bmi, pardm, insoroha, hba1c, C, A, M, Gene) 

dataset_type2 <- full_join(MYDIABETES_type2, 
                           UNITED_type2p, 
                           by = c("study",
                                  "agerec", 
                                  "agedx", 
                                  "sex", 
                                  "bmi", 
                                  "pardm", 
                                  "insoroha", 
                                  "hba1c", 
                                  "C", 
                                  "A", 
                                  "M", 
                                  "Gene")) %>%
  mutate(
    Gene = ifelse(Gene == "", NA, Gene), 
    biomark_status = ifelse(is.na(biomark_status), 
                            ifelse(C == 1 & A == 0, 0, 1),
                            biomark_status)) 

external <- full_join(dataset_type1, dataset_type2)

cc_type1 %>%
  summarise(mean = mean(duration),
            median = median(duration))
CC %>%
  summarise(mean = mean(duration),
            median = median(duration))
CC %>%
  summarise(
    n = sum(duration < 5), 
    perc = n/n()
  )

CC %>%
  summarise(
    n = sum(duration < 3), 
    perc = n/n(),
    min = min(duration)
  )

CC %>%
  filter(duration < 0) %>%
  select(agedx, agerec, typedm, duration)
cc_type2 %>%
  summarise(mean = mean(duration),
            median = median(duration))
ggplot(cc_type1, aes(x= duration)) +
  geom_histogram()
ggplot(cc_type1, aes(x= duration)) +
  geom_boxplot()
ggplot(cc_type2, aes(x= duration)) +
  geom_boxplot()

UNITED_type1 %>%
  summarise(mean = mean(durationfinal),
            median = median(durationfinal))
UNITED %>%
  summarise(mean = mean(durationfinal),
            median = median(durationfinal))

UNITED  %>%
  filter(durationfinal < 0) %>%
  select(agedx, agerec)

UNITED %>%
  summarise(
    n = sum(durationfinal < 5), 
    perc = n/n()
  )
UNITED %>%
  summarise(
    n = sum(durationfinal < 3), 
    perc = n/n(),
    min = min(durationfinal)
  )
UNITED_type2 %>%
  summarise(mean = mean(durationfinal),
            median = median(durationfinal))
ggplot(UNITED_type1, aes(x= durationfinal)) +
  geom_boxplot()
ggplot(UNITED_type2, aes(x= durationfinal)) +
  geom_boxplot()
ggplot(UNITED, aes(x= durationfinal)) +
  geom_boxplot()

ggplot(UNITED, aes(x= durationfinal)) +
  geom_boxplot()


UNITED_duration <- UNITED %>%
  mutate(
    perc_5years = (sum(durationfinal < 5)/n())*100, 
    perc_3years = (sum(durationfinal < 3)/n())*100
  ) %>%
  ggplot(aes(x= durationfinal)) +
  geom_boxplot() + 
  geom_vline(xintercept = 5) +  
  geom_vline(xintercept = 3) +
  theme_bw() +
  xlab("Duration of diabetes at time of study (years)")

CC_duration <- CC %>%
  mutate(
    perc_5years = (sum(duration < 5)/n())*100, 
    perc_3years = (sum(duration < 3)/n())*100
  ) %>%
  ggplot(aes(x= duration)) +
  geom_boxplot() + 
  geom_vline(xintercept = 5) +  
  geom_vline(xintercept = 3) +
  theme_bw() +
  xlab("Duration of diabetes at time of study (years)")

external %>%
  mutate(
    durationfinal = ifelse(is.na(durationfinal), 
                           agerec - agedx, 
                           durationfinal)) %>%
  summarise(
    perc_5years = (sum(durationfinal < 5)/n())*100, 
    perc_3years = (sum(durationfinal < 3)/n())*100, 
    mean = mean(durationfinal),
    median = median(durationfinal),
    min = min(durationfinal))

external  %>%
  mutate(
    durationfinal = ifelse(is.na(durationfinal), 
                           agerec - agedx, 
                           durationfinal)) %>%
  filter(durationfinal < 0) %>%
  select(agedx, agerec, durationfinal, study, MY_ID)

external_duration <- external %>%
  mutate(
    durationfinal = ifelse(is.na(durationfinal), 
                           agerec - agedx, 
                           durationfinal),
    perc_5years = (sum(durationfinal < 5)/n())*100, 
    perc_3years = (sum(durationfinal < 3)/n())*100
  ) %>%
  ggplot(aes(x= durationfinal)) +
  geom_boxplot() + 
  geom_vline(xintercept = 5) +  
  geom_vline(xintercept = 3) +
  theme_bw() +
  xlab("Duration of diabetes at time of study (years)")

data_duration_plot <- patchwork::wrap_plots(
  CC_duration, 
  UNITED_duration,
  external_duration,
  ncol = 1, nrow = 3
) + patchwork::plot_annotation(tag_levels = list(c("A", "B", "C"))) &
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 11), 
    plot.margin = margin(10,10,10,10)
  ) 

pdf("Data_duration.pdf", height = 12, width = 18)
print(data_duration_plot)
dev.off()


