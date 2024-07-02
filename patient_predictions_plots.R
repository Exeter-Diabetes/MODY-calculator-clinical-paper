#:--------------------------------------------------------
#
# In this file we plot the probabilities for individuals
#
#:--------------------------------------------------------

# load libraries
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggforce)
library(pROC)

# load files required
predictions_dataset.UNITED_type1_with_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_with_T.rds")
predictions_dataset.referral_type1_with_T <- readRDS("model_predictions/predictions_dataset.referral_type1_with_T.rds")
predictions_dataset.UNITED_type2_new <- readRDS("model_predictions/predictions_dataset.UNITED_type2_new.rds")

probabilities_under_30 <- read.delim("data/mody_probabilities_under_30s.txt")

# load functions needed
source("data/create_data.R")
calculate_thresholds_diagnostics <- function(response, prediction, unique = FALSE) {
  
  ### threshold number to test
  thresholds <- seq(0, 1, 0.0001)
  
  matrix_thresholds <- matrix(0, nrow = length(thresholds), ncol = 8)
  matrix_thresholds[,1] <- thresholds
  
  preds1 <- data.frame(
    prediction = prediction,
    response = response
  )
  
  for (i in 1:length(thresholds)) {
    # original dataset
    original_data <- preds1
    
    # checking which patients are above the threshold
    interim_above <- preds1 %>%
      mutate(threshold_up = ifelse(prediction > thresholds[i], 1, 0)) %>%
      filter(threshold_up == 1)
    
    interim_below <- preds1 %>%
      mutate(threshold_up = ifelse(prediction > thresholds[i], 1, 0)) %>%
      filter(threshold_up == 0)
    
    # calculating the sensitivity of the MODY cases
    matrix_thresholds[i, 2] <- length(which(interim_above$response == 1))/length(which(original_data$response == 1)) * 100
    
    # calculating number to test
    if (length(which(interim_above$response == 1)) > 0) {
      matrix_thresholds[i, 3] <- ceiling(nrow(interim_above)/length(which(interim_above$response == 1)))
    } else {matrix_thresholds[i, 3] <- as.numeric(NA)}
    
    # Pick-up rate
    matrix_thresholds[i, 4] <- length(which(interim_above$response == 1))/nrow(interim_above) * 100
    
    # calculating specificity
    matrix_thresholds[i, 5] <- length(which(interim_below$response == 0))/length(which(original_data$response == 0)) * 100
    
    # Perc of patients tested
    matrix_thresholds[i, 6] <- (nrow(interim_above)/nrow(original_data))*100
    
    # number of MODY cases missed
    matrix_thresholds[i, 7] <- interim_below %>%
      select(response) %>%
      sum(na.rm = TRUE) %>%
      unlist()
    
    # number of patients tested
    matrix_thresholds[i, 8] <- interim_above %>%
      nrow() %>%
      unlist()
    
  }
  
  ## making the matrix a data.frame()
  matrix_thresholds <- matrix_thresholds %>% 
    as.data.frame() %>%
    set_names(c("Thresholds", "Sensitivity", "NTT", "Pick-up rate", "Specificity", "N-Tested", "Cases Missed", "Patients Tested")) %>%
    arrange(desc(Sensitivity), desc(NTT), desc(`Pick-up rate`))
  
  if (unique == TRUE) {
    ## select unique combinations of sensitivity and NTT (only the first occurance)
    matrix_thresholds <- matrix_thresholds %>%
      slice(which(duplicated(matrix_thresholds %>% select(-Thresholds, -`Pick-up rate`)) == FALSE))
  }
  
  return(matrix_thresholds)
  
}



## Load population representative dataset
dataset.UNITED_type1 <- create_data(dataset = "united t1d") %>%
  
  ## if MODY testing missing, change to 0
  mutate(M = ifelse(is.na(M), 0, M))

dataset.UNITED_type2 <- create_data(dataset = "united t2d") %>%
  
  ## if MODY testing missing, change to 0
  mutate(M = ifelse(is.na(M), 0, M))


## Load referrals dataset
### load in Referral repository
current_path <- getwd()
setwd(dir = paste0(current_path, "/data/"))
#### download repository to get the latest updates
download.file(
  url = "https://github.com/Exeter-Diabetes/Julieanne-Pedro-MODY-Referrals/archive/master.zip",
  destfile = "meetingsR-master.zip"
)
#### unzip github repository
unzip(zipfile = "meetingsR-master.zip")
#### remove zipped file
file.remove("meetingsR-master.zip")
#### reverse working directory
setwd(current_path)
#### forget unnecessary variables
rm(current_path, create_data)

### load datasets
library(readxl)
source("data/Julieanne-Pedro-MODY-Referrals-main/data_formatting.R")

dataset.referral_type1 <- formatting(
  dataset = read_excel("data/mody_35yrs_08_02_2024.xlsx", guess_max = 1000000),
  read_csv("data/mmoct11.csv"),
  ethnicity_groups = read_excel("data/Julieanne-Pedro-MODY-Referrals-main/ethnicity_groups_ngs.xlsx", na = "null"),
  ethnicity_groups_genetics = read_excel("data/Julieanne-Pedro-MODY-Referrals-main/ancestry_ngs_kp.xlsx", na = "null"),
  ethnicity_labels_stated = read_excel("data/Julieanne-Pedro-MODY-Referrals-main/ethnicity_labels_stated.xlsx", na = "null"),
  ethnicity_labels_genetics = read_excel("data/Julieanne-Pedro-MODY-Referrals-main/ethnicity_labels_genetics.xlsx", na = "null"),
  biomarkers = "any",
  diagnosis = FALSE,
  type = "Type 1",
  ethnicity = "White",
  proband = "Proband",
  gene_variable = FALSE
) %>%
  drop_na(c("sex", "bmi", "agedx", "hba1c", "pardm", "agerec")) %>%   # what should be done about missing MODY testing?
  mutate(T = ifelse(C == 0 | A == 1, 1, 0)) # T is 1 if Cn or Ap

#### remove downloaded folder
unlink("data/Julieanne-Pedro-MODY-Referrals-main", recursive = TRUE)


#:--------------------------------------------

# Table information

## Insulin treated
thresholds_UNITED_t1d <- calculate_thresholds_diagnostics(dataset.UNITED_type1$M, predictions_dataset.UNITED_type1_with_T$prob)

### 5%
thresholds_UNITED_t1d %>%
  filter(Thresholds == 0.05) %>% arrange(Thresholds) %>% head()

probabilities_under_30 %>%
  filter(which_equation == "t1" & pedro_prob > 5) %>%
  nrow()

### 10%
thresholds_UNITED_t1d %>%
  filter(Thresholds == 0.1) %>% arrange(Thresholds)  %>% head()

probabilities_under_30 %>%
  filter(which_equation == "t1" & pedro_prob > 10) %>%
  nrow()

### 20%
thresholds_UNITED_t1d %>%
  filter(Thresholds == 0.2) %>% arrange(Thresholds)  %>% head()

probabilities_under_30 %>%
  filter(which_equation == "t1" & pedro_prob > 20) %>%
  nrow()

### 30%
thresholds_UNITED_t1d %>%
  filter(Thresholds == 0.3) %>% arrange(Thresholds)  %>% head()

probabilities_under_30 %>%
  filter(which_equation == "t1" & pedro_prob > 30) %>%
  nrow()


## Non-insulin treated
thresholds_UNITED_t2d <- calculate_thresholds_diagnostics(dataset.UNITED_type2$M, predictions_dataset.UNITED_type2_new$prob)

### 5%
thresholds_UNITED_t2d %>%
  filter(Thresholds == 0.05) %>% arrange(Thresholds)  %>% head()

probabilities_under_30 %>%
  filter(which_equation == "t2" & pedro_prob > 5) %>%
  nrow()

### 10%
thresholds_UNITED_t2d %>%
  filter(Thresholds == 0.1) %>% arrange(Thresholds)  %>% head()

probabilities_under_30 %>%
  filter(which_equation == "t2" & pedro_prob > 10) %>%
  nrow()

### 20%
thresholds_UNITED_t2d %>%
  filter(Thresholds == 0.2) %>% arrange(Thresholds)  %>% head()

probabilities_under_30 %>%
  filter(which_equation == "t2" & pedro_prob > 20) %>%
  nrow()

### 25%
thresholds_UNITED_t2d %>%
  filter(Thresholds == 0.25) %>% arrange(Thresholds)  %>% head()

probabilities_under_30 %>%
  filter(which_equation == "t2" & pedro_prob > 25) %>%
  nrow()

### 30%
thresholds_UNITED_t2d %>%
  filter(Thresholds == 0.3) %>% arrange(Thresholds)  %>% head()

probabilities_under_30 %>%
  filter(which_equation == "t2" & pedro_prob > 30) %>%
  nrow()


#:--------------------------------------------


##: All patients in UNITED

thresholds_UNITED_t1d <- calculate_thresholds_diagnostics(dataset.UNITED_type1$M, predictions_dataset.UNITED_type1_with_T$prob)


# combine all the necessary columns
plot_probabilities_united <- cbind(
  mody = dataset.UNITED_type1$M,
  prob = predictions_dataset.UNITED_type1_with_T$prob,
  LCI = predictions_dataset.UNITED_type1_with_T$LCI,
  UCI = predictions_dataset.UNITED_type1_with_T$UCI
) %>%
  
  # turn columns into a data.frame
  as.data.frame() %>%
  
  # change mody column to factor
  mutate(
    mody = factor(mody, levels = c(1, 0), labels = c("MODY", "Non-MODY"))
  ) %>%
  
  # sort patients based on the mean prob (per MODY status)
  arrange(desc(prob)) %>%
  
  # add a ranking variable
  mutate(rank = 1:n()) %>%
  
  # plot probabilities
  ggplot() +
  geom_hline(yintercept = 0, colour = "black") +
  geom_pointrange(aes(x = rank, y = prob, ymin = LCI, ymax = UCI), colour = "grey") +
  geom_point(aes(x = rank, y = prob)) +
  geom_hline(yintercept = 0.1468, colour = "#E69F00") + ## optimal cut-off
  geom_hline(yintercept = 0.1, colour = "#D55E00") +
  geom_hline(yintercept = 0.25, colour = "#56B4E9") +
  scale_y_continuous("MODY Probability", labels = scales::percent, breaks = seq(0, 1, 0.1)) +
  scale_x_continuous("Ranked patients") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18)
  )



cbind(
  mody = dataset.UNITED_type1$M,
  prob = predictions_dataset.UNITED_type1_with_T$prob
) %>%
  as.data.frame() %>%
  filter(prob > 0.25) %>%
  # filter(prob > 0.1468) %>%
  select(mody) %>%
  unlist() %>%
  table()


# combine all the necessary columns
plot_probabilities_united_ppv <- cbind(
  Threshold = c(
    "25%\nn = 6 MODY = 1", "25%\nn = 6 MODY = 1", 
    "Optimal: 14.7%\nn = 15 MODY = 5", "Optimal: 14.7%\nn = 15 MODY = 5",
    "10%\nn = 23 MODY = 6", "10%\nn = 23 MODY = 6",
    "0%\nn = 1171 MODY = 7", "0%\nn = 1171 MODY = 7"
  ),
  mody = c("Non-MODY", "MODY", "Non-MODY", "MODY", "Non-MODY", "MODY", "Non-MODY", "MODY"),
  count = c(
    5, 1,
    10, 5,
    17, 6,
    1164, 7 
  )
) %>%
  as.data.frame() %>%
  mutate(
    Threshold = factor(Threshold, levels = c("0%\nn = 1171 MODY = 7", "10%\nn = 23 MODY = 6", "Optimal: 14.7%\nn = 15 MODY = 5", "25%\nn = 6 MODY = 1")),
    count = as.numeric(count),
    mody = factor(mody, levels = c("Non-MODY", "MODY"))
  ) %>%
  ggplot() +
  geom_bar(
    aes(fill = mody, y = Threshold, x = count), colour = "black", position = "fill", stat = "identity"
  ) +
  scale_fill_manual(values = c("black", "white"), breaks=c("MODY", "Non-MODY")) +
  scale_x_continuous(labels = scales::percent, "Positive Predictive Value", breaks = seq(0, 1, 0.1)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title=element_blank(),
    axis.title = element_text(size = 18),
    axis.text.y = element_text(size = 15, colour = c("black", "#D55E00", "#E69F00", "#56B4E9")),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 18)
  )




plot_probabilities_combined_united <- patchwork::wrap_plots(
  
  plot_probabilities_united,
  plot_probabilities_united_ppv,
  ncol = 1
  
) + patchwork::plot_annotation(tag_levels = 'A')




##: All patients in referrals

thresholds <- calculate_thresholds_diagnostics(dataset.referral_type1$M, predictions_dataset.referral_type1_with_T$prob)


# combine all the necessary columns
plot_probabilities_referral <- cbind(
  mody = dataset.referral_type1$M,
  prob = predictions_dataset.referral_type1_with_T$prob,
  LCI = predictions_dataset.referral_type1_with_T$LCI,
  UCI = predictions_dataset.referral_type1_with_T$UCI
) %>%
  
  # turn columns into a data.frame
  as.data.frame() %>%
  
  # change mody column to factor
  mutate(
    mody = factor(mody, levels = c(1, 0), labels = c("MODY", "Non-MODY"))
  ) %>%
  
  # sort patients based on the mean prob (per MODY status)
  arrange(desc(prob)) %>%
  
  # add a ranking variable
  mutate(rank = 1:n()) %>%
  
  # plot probabilities
  ggplot() +
  geom_hline(yintercept = 0, colour = "black") +
  geom_pointrange(aes(x = rank, y = prob, ymin = LCI, ymax = UCI), colour = "grey") +
  geom_point(aes(x = rank, y = prob)) +
  geom_hline(yintercept = 0.63, colour = "#E69F00") + ## optimal cut-off referral
  geom_hline(yintercept = 0.1468, colour = "#CC79A7") +  ## optimal cut-off UNITED
  geom_hline(yintercept = 0.1, colour = "#D55E00") +
  geom_hline(yintercept = 0.25, colour = "#56B4E9") +
  geom_hline(yintercept = 0.5, colour = "#009E73") +
  geom_hline(yintercept = 0.75, colour = "#0072B2") +
  scale_y_continuous("MODY Probability", labels = scales::percent, breaks = seq(0, 1, 0.1)) +
  scale_x_continuous("Ranked patients") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18)
  )



# cbind(
#   mody = dataset.referral_type1$M,
#   prob = predictions_dataset.referral_type1_with_T$prob
# ) %>%
#   as.data.frame() %>%
#   # filter(prob > 0.1) %>%
#   filter(prob > 0.1468) %>%
#   select(mody) %>%
#   unlist() %>%
#   table()





# combine all the necessary columns
plot_probabilities_referral_ppv <- cbind(
  Threshold = c(
    "75%\nn = 4 MODY = 1", "75%\nn = 4 MODY = 1",
    "Ref Optimal: 63%\nn = 45 MODY = 13", "Ref Optimal: 63%\nn = 45 MODY = 13",
    "50%\nn = 81 MODY = 21", "50%\nn = 81 MODY = 21", 
    "25%\nn = 213 MODY = 56", "25%\nn = 213 MODY = 56", 
    "UNITED Optimal: 14.7%\nn = 384 MODY = 96", "UNITED Optimal: 14.7%\nn = 384 MODY = 96",
    "10%\nn = 492 MODY = 118", "10%\nn = 492 MODY = 118",
    "0%\nn = 1754 MODY = 229", "0%\nn = 1754 MODY = 229"
  ),
  mody = c("Non-MODY", "MODY", "Non-MODY", "MODY", "Non-MODY", "MODY", "Non-MODY", "MODY", "Non-MODY", "MODY", "Non-MODY", "MODY", "Non-MODY", "MODY"),
  count = c(
    3, 1,
    32, 13,
    63, 21,
    157, 56,
    288, 96,
    374, 118,
    1525, 229
  )
) %>%
  as.data.frame() %>%
  mutate(
    Threshold = factor(Threshold, levels = c("0%\nn = 1754 MODY = 229", "10%\nn = 492 MODY = 118", "UNITED Optimal: 14.7%\nn = 384 MODY = 96", "25%\nn = 213 MODY = 56", "50%\nn = 81 MODY = 21", "Ref Optimal: 63%\nn = 45 MODY = 13", "75%\nn = 4 MODY = 1")),
    count = as.numeric(count),
    mody = factor(mody, levels = c("Non-MODY", "MODY"))
  ) %>%
  ggplot() +
  geom_bar(
    aes(fill = mody, y = Threshold, x = count), colour = "black", position = "fill", stat = "identity"
  ) +
  scale_fill_manual(values = c("black", "white"), breaks=c("MODY", "Non-MODY")) +
  scale_x_continuous(labels = scales::percent, "Positive Predictive Value", breaks = seq(0, 1, 0.1)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title=element_blank(),
    axis.title = element_text(size = 18),
    axis.text.y = element_text(size = 15, colour = c("black", "#D55E00", "#CC79A7", "#56B4E9", "#009E73", "#E69F00", "#0072B2")),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 18)
  )




plot_probabilities_combined_referral <- patchwork::wrap_plots(
  
  plot_probabilities_referral,
  plot_probabilities_referral_ppv,
  ncol = 1
  
) + patchwork::plot_annotation(tag_levels = 'A')




#:--------------------------------------------------------

# Referrals with biomarker info


# calculate PPV for thresholds in those with biomarker info
biomarker_available <- dataset.referral_type1 %>%
  cbind(prob = predictions_dataset.referral_type1_with_T$prob) %>%
  filter(!is.na(T))

thresholds <- calculate_thresholds_diagnostics(biomarker_available$M, biomarker_available$prob)


# combine all the necessary columns
plot_probabilities_referral_biomarker_info <- data.frame(mody = dataset.referral_type1$M) %>%
  
  cbind(
    T = dataset.referral_type1$T,
    prob = predictions_dataset.referral_type1_with_T$prob,
    LCI = predictions_dataset.referral_type1_with_T$LCI,
    UCI = predictions_dataset.referral_type1_with_T$UCI
  ) %>%
  
  # keep only those with biomarker info
  filter(!is.na(T)) %>%
  
  # turn columns into a data.frame
  as.data.frame() %>%
  
  # change mody column to factor
  mutate(
    mody = factor(mody, levels = c(1, 0), labels = c("MODY", "Non-MODY"))
  ) %>%
  
  # sort patients based on the mean prob (per MODY status)
  arrange(desc(prob)) %>%
  
  # add a ranking variable
  mutate(rank = 1:n()) %>%
  
  # plot probabilities
  ggplot() +
  geom_hline(yintercept = 0, colour = "black") +
  geom_pointrange(aes(x = rank, y = prob, ymin = LCI, ymax = UCI), colour = "grey") +
  geom_point(aes(x = rank, y = prob)) +
  geom_hline(yintercept = 0.63, colour = "#E69F00") + ## optimal cut-off
  geom_hline(yintercept = 0.1468, colour = "#CC79A7") +  ## optimal cut-off UNITED
  geom_hline(yintercept = 0.1, colour = "#D55E00") +
  geom_hline(yintercept = 0.25, colour = "#56B4E9") +
  geom_hline(yintercept = 0.5, colour = "#009E73") +
  geom_hline(yintercept = 0.75, colour = "#0072B2") +
  scale_y_continuous("MODY Probability", labels = scales::percent, breaks = seq(0, 1, 0.1)) +
  scale_x_continuous("Ranked patients") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18)
  )



# dataset.referral_type1 %>%
#   cbind(prob = predictions_dataset.referral_type1_with_T$prob) %>%
#   filter(!is.na(T)) %>%
#   as.data.frame() %>%
#   # filter(prob > 0.10) %>%
#   filter(prob > 0.1468) %>%
#   select(M) %>%
#   unlist() %>%
#   table()





# combine all the necessary columns
plot_probabilities_referral_ppv_biomarker_info <- cbind(
  Threshold = c(
    "75%\nn = 4 MODY = 1", "75%\nn = 4 MODY = 1",
    "Ref Optimal: 63%\nn = 44 MODY = 13", "Ref Optimal: 63%\nn = 44 MODY = 13",
    "50%\nn = 84 MODY = 21", "50%\nn = 84 MODY = 21", 
    "25%\nn = 192 MODY = 49", "25%\nn = 192 MODY = 49", 
    "UNITED Optimal: 14.7%\nn = 264 MODY = 61", "UNITED Optimal: 14.7%\nn = 264 MODY = 61",
    "10%\nn = 307 MODY = 69", "10%\nn = 307 MODY = 69", 
    "0%\nn = 992 MODY = 111", "0%\nn = 992 MODY = 111"
  ),
  mody = c("Non-MODY", "MODY", "Non-MODY", "MODY", "Non-MODY", "MODY", "Non-MODY", "MODY", "Non-MODY", "MODY", "Non-MODY", "MODY", "Non-MODY", "MODY"),
  count = c(
    3, 1,
    31, 13,
    63, 21,
    143, 49,
    203, 61,
    238, 69,
    881, 111
  )
) %>%
  as.data.frame() %>%
  mutate(
    Threshold = factor(Threshold, levels = c("0%\nn = 992 MODY = 111", "10%\nn = 307 MODY = 69", "UNITED Optimal: 14.7%\nn = 264 MODY = 61", "25%\nn = 192 MODY = 49", "50%\nn = 84 MODY = 21", "Ref Optimal: 63%\nn = 44 MODY = 13", "75%\nn = 4 MODY = 1")),
    count = as.numeric(count),
    mody = factor(mody, levels = c("Non-MODY", "MODY"))
  ) %>%
  ggplot() +
  geom_bar(
    aes(fill = mody, y = Threshold, x = count), colour = "black", position = "fill", stat = "identity"
  ) +
  scale_fill_manual(values = c("black", "white"), breaks=c("MODY", "Non-MODY")) +
  scale_x_continuous(labels = scales::percent, "Positive Predictive Value", breaks = seq(0, 1, 0.1)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title=element_blank(),
    axis.title = element_text(size = 18),
    axis.text.y = element_text(size = 15, colour = c("black", "#D55E00", "#CC79A7", "#56B4E9", "#009E73", "#E69F00", "#0072B2")),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 18)
  )




plot_probabilities_combined_biomarker_info <- patchwork::wrap_plots(
  
  plot_probabilities_referral_biomarker_info,
  plot_probabilities_referral_ppv_biomarker_info,
  ncol = 1
  
) + patchwork::plot_annotation(tag_levels = 'A')












