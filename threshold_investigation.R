#:--------------------------------------------------------
#
# In this file we make waffle plots for threshold investigation
#
#:--------------------------------------------------------

# load libraries
library(tidyverse)
library(ggplot2)
library(waffle)
library(patchwork)


# load files required
predictions_dataset.UNITED_type1_old <- readRDS("model_predictions/predictions_dataset.UNITED_type1_old.rds")
predictions_dataset.UNITED_type1_with_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_with_T.rds")
predictions_dataset.UNITED_type1_no_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_no_T.rds")
predictions_dataset.UNITED_type2_new <- readRDS("model_predictions/predictions_dataset.UNITED_type2_new.rds")

predictions_dataset.referral_type1_old <- readRDS("model_predictions/predictions_dataset.referral_type1_old.rds")
predictions_dataset.referral_type1_with_T <- readRDS("model_predictions/predictions_dataset.referral_type1_with_T.rds")
predictions_dataset.referral_type1_no_T <- readRDS("model_predictions/predictions_dataset.referral_type1_no_T.rds")
predictions_dataset.referral_type2_new <- readRDS("model_predictions/predictions_dataset.referral_type2_new.rds")


# load functions needed
source("data/create_data.R")


## Load population representative dataset
dataset.UNITED_type1 <- create_data(dataset = "united t1d") %>%
  
  ## if MODY testing missing, change to 0
  mutate(M = ifelse(is.na(M), 0, M))

dataset.UNITED_type2 <- create_data(dataset = "united t2d")



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

dataset.referral_type2 <- formatting(
  dataset = read_excel("data/mody_35yrs_08_02_2024.xlsx", guess_max = 1000000), 
  read_csv("data/mmoct11.csv"), 
  ethnicity_groups = read_excel("data/Julieanne-Pedro-MODY-Referrals-main/ethnicity_groups_ngs.xlsx", na = "null"), 
  ethnicity_groups_genetics = read_excel("data/Julieanne-Pedro-MODY-Referrals-main/ancestry_ngs_kp.xlsx", na = "null"), 
  ethnicity_labels_stated = read_excel("data/Julieanne-Pedro-MODY-Referrals-main/ethnicity_labels_stated.xlsx", na = "null"), 
  ethnicity_labels_genetics = read_excel("data/Julieanne-Pedro-MODY-Referrals-main/ethnicity_labels_genetics.xlsx", na = "null"), 
  diagnosis = FALSE, 
  type = "Type 2", 
  ethnicity = "White", 
  proband = "Proband", 
  gene_variable = FALSE
) %>%
  drop_na(c("sex", "bmi", "agedx", "insoroha", "hba1c", "pardm", "agerec"))    # what should be done about missing MODY testing?

#### remove downloaded folder
unlink("data/Julieanne-Pedro-MODY-Referrals-main", recursive = TRUE)



#:--------------------------------------------

plot_UNITED_type1_no_T <- patchwork::wrap_plots(
  
  dataset.UNITED_type1 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.UNITED_type1_no_T$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.1, "Over", "Under")) %>%
    select(threshold) %>%
    group_by(threshold) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(threshold = ifelse(threshold == "Under", paste0("<10%\nn=", n), paste0(">10%\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = threshold, values = n), colour = "white") +
    theme_void() +
    ggtitle(paste0("Insulin treated patients\nn=", nrow(predictions_dataset.UNITED_type1_no_T))) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ),
  
  dataset.UNITED_type1 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.UNITED_type1_no_T$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.1, "Over", "Under")) %>%
    filter(threshold == "Over") %>%
    select(M) %>%
    group_by(M) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    rbind(cbind(M = 1, n = 0)) %>%
    mutate(M = ifelse(M == 0, paste0("Non-MODY\nn=", n), paste0("MODY\nn=", n))) %>%   #, limits = c("Non-MODY\nn=8", "MODY\nn=0")
    ggplot() +
    geom_waffle(aes(fill = M, values = n), colour = "white", show.legend = TRUE) +
    scale_fill_manual(values = c("Non-MODY\nn=8" = "deepskyblue4", "MODY\nn=0" = "#00BFC4"), limits = c("Non-MODY\nn=8", "MODY\nn=0")) +
    theme_void() +
    ggtitle(paste0("Insulin treated patients\nwith MODY probability >10%")) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ),
  
  dataset.UNITED_type1 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.UNITED_type1_no_T$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.1, "Over", "Under")) %>%
    filter(threshold == "Under") %>%
    select(M) %>%
    group_by(M) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(M = ifelse(M == 0, paste0("Non-MODY\nn=", n), paste0("MODY\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = M, values = n), colour = "white") +
    scale_fill_manual(values = c("red", "red4")) +
    theme_void() +
    ggtitle(paste0("Insulin treated patients\nwith MODY probability <10%")) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ),
  
  nrow = 1
  
) + 
  patchwork::plot_annotation(
    title = "UNITED: insulin - clinical features (10% threshold)"
  ) &
  theme(
    legend.position = "bottom"
  )


plot_UNITED_type1_with_T <- patchwork::wrap_plots(
  
  dataset.UNITED_type1 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.UNITED_type1_with_T$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.1, "Over", "Under")) %>%
    select(threshold) %>%
    group_by(threshold) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(threshold = ifelse(threshold == "Under", paste0("<10%\nn=", n), paste0(">10%\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = threshold, values = n), colour = "white") +
    theme_void() +
    ggtitle(paste0("Insulin treated patients\nn=", nrow(predictions_dataset.UNITED_type1_with_T))) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ),
  
  dataset.UNITED_type1 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.UNITED_type1_with_T$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.1, "Over", "Under")) %>%
    filter(threshold == "Over") %>%
    select(M) %>%
    group_by(M) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(M = ifelse(M == 0, paste0("Non-MODY\nn=", n), paste0("MODY\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = M, values = n), colour = "white", show.legend = TRUE) +
    scale_fill_manual(values = c("#00BFC4", "deepskyblue4")) +
    theme_void() +
    ggtitle(paste0("Insulin treated patients\nwith MODY probability >10%")) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ),
  
  dataset.UNITED_type1 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.UNITED_type1_with_T$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.1, "Over", "Under")) %>%
    filter(threshold == "Under") %>%
    select(M) %>%
    group_by(M) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(M = ifelse(M == 0, paste0("Non-MODY\nn=", n), paste0("MODY\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = M, values = n), colour = "white") +
    scale_fill_manual(values = c("red", "red4")) +
    theme_void() +
    ggtitle(paste0("Insulin treated patients\nwith MODY probability <10%")) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ),
  
  nrow = 1
  
) +
  patchwork::plot_annotation(
    title = "UNITED: insulin - biomarkers (10% threshold)"
  ) &
  theme(
    legend.position = "bottom"
  )


plot_UNITED_type1_old <- patchwork::wrap_plots(
  
  dataset.UNITED_type1 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.UNITED_type1_old$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.1, "Over", "Under")) %>%
    select(threshold) %>%
    group_by(threshold) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(threshold = ifelse(threshold == "Under", paste0("<10%\nn=", n), paste0(">10%\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = threshold, values = n), colour = "white") +
    theme_void() +
    ggtitle(paste0("Insulin treated patients\nn=", nrow(predictions_dataset.UNITED_type1_old))) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ),
  
  dataset.UNITED_type1 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.UNITED_type1_old$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.1, "Over", "Under")) %>%
    filter(threshold == "Over") %>%
    select(M) %>%
    group_by(M) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(M = ifelse(M == 0, paste0("Non-MODY\nn=", n), paste0("MODY\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = M, values = n), colour = "white", show.legend = TRUE) +
    scale_fill_manual(values = c("#00BFC4", "deepskyblue4")) +
    theme_void() +
    ggtitle(paste0("Insulin treated patients\nwith MODY probability >10%")) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ),
  
  dataset.UNITED_type1 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.UNITED_type1_old$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.1, "Over", "Under")) %>%
    filter(threshold == "Under") %>%
    select(M) %>%
    group_by(M) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(M = ifelse(M == 0, paste0("Non-MODY\nn=", n), paste0("MODY\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = M, values = n), colour = "white") +
    scale_fill_manual(values = c("red", "red4")) +
    theme_void() +
    ggtitle(paste0("Insulin treated patients\nwith MODY probability <10%")) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ),
  
  nrow = 1
  
) +
  patchwork::plot_annotation(
    title = "UNITED: insulin - old (10% threshold)"
  ) &
  theme(
    legend.position = "bottom"
  )



plot_UNITED_type2 <- patchwork::wrap_plots(
  
  dataset.UNITED_type2 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.UNITED_type2_new$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.25, "Over", "Under")) %>%
    select(threshold) %>%
    group_by(threshold) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(threshold = ifelse(threshold == "Under", paste0("<25%\nn=", n), paste0(">25%\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = threshold, values = n), colour = "white") +
    theme_void() +
    ggtitle(paste0("Non-insulin treated patients\nn=", nrow(predictions_dataset.UNITED_type2_new))) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ),
  
  dataset.UNITED_type2 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.UNITED_type2_new$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.25, "Over", "Under")) %>%
    filter(threshold == "Over") %>%
    select(M) %>%
    group_by(M) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(M = ifelse(M == 0, paste0("Non-MODY\nn=", n), paste0("MODY\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = M, values = n), colour = "white", show.legend = TRUE) +
    scale_fill_manual(values = c("#00BFC4", "deepskyblue4")) +
    theme_void() +
    ggtitle(paste0("Non-insulin treated patients\nwith MODY probability >25%")) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ),
  
  dataset.UNITED_type2 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.UNITED_type2_new$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.25, "Over", "Under")) %>%
    filter(threshold == "Under") %>%
    select(M) %>%
    group_by(M) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(M = ifelse(M == 0, paste0("Non-MODY\nn=", n), paste0("MODY\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = M, values = n), colour = "white") +
    scale_fill_manual(values = c("red", "red4")) +
    theme_void() +
    ggtitle(paste0("Non-insulin treated patients\nwith MODY probability <25%")) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ),
  
  nrow = 1
  
) +
  patchwork::plot_annotation(
    title = "UNITED: non-insulin - clinical features (25% threshold)"
  ) &
  theme(
    legend.position = "bottom"
  )


plot_referral_type1_no_T <- patchwork::wrap_plots(
  
  dataset.referral_type1 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.referral_type1_no_T$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.1, "Over", "Under")) %>%
    select(threshold) %>%
    group_by(threshold) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(threshold = ifelse(threshold == "Under", paste0("<10%\nn=", n), paste0(">10%\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = threshold, values = n), colour = "white") +
    theme_void() +
    ggtitle(paste0("Insulin treated patients\nn=", nrow(predictions_dataset.referral_type1_no_T))) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ),
  
  dataset.referral_type1 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.referral_type1_no_T$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.1, "Over", "Under")) %>%
    filter(threshold == "Over") %>%
    select(M) %>%
    group_by(M) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(M = ifelse(M == 0, paste0("Non-MODY\nn=", n), paste0("MODY\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = M, values = n), colour = "white", show.legend = TRUE) +
    scale_fill_manual(values = c("#00BFC4", "deepskyblue4")) +
    theme_void() +
    ggtitle(paste0("Insulin treated patients\nwith MODY probability >10%")) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ),
  
  dataset.referral_type1 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.referral_type1_no_T$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.1, "Over", "Under")) %>%
    filter(threshold == "Under") %>%
    select(M) %>%
    group_by(M) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(M = ifelse(M == 0, paste0("Non-MODY\nn=", n), paste0("MODY\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = M, values = n), colour = "white") +
    scale_fill_manual(values = c("red", "red4")) +
    theme_void() +
    ggtitle(paste0("Insulin treated patients\nwith MODY probability <10%")) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ),
  
  nrow = 1
  
) +
  patchwork::plot_annotation(
    title = "Referrals: insulin - clinical features (10% threshold)"
  ) &
  theme(
    legend.position = "bottom"
  )


plot_referral_type1_with_T <- patchwork::wrap_plots(
  
  dataset.referral_type1 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.referral_type1_with_T$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.1, "Over", "Under")) %>%
    select(threshold) %>%
    group_by(threshold) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(threshold = ifelse(threshold == "Under", paste0("<10%\nn=", n), paste0(">10%\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = threshold, values = n), colour = "white") +
    theme_void() +
    ggtitle(paste0("Insulin treated patients\nn=", nrow(predictions_dataset.referral_type1_with_T))) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ),
  
  dataset.referral_type1 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.referral_type1_with_T$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.1, "Over", "Under")) %>%
    filter(threshold == "Over") %>%
    select(M) %>%
    group_by(M) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(M = ifelse(M == 0, paste0("Non-MODY\nn=", n), paste0("MODY\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = M, values = n), colour = "white", show.legend = TRUE) +
    scale_fill_manual(values = c("#00BFC4", "deepskyblue4")) +
    theme_void() +
    ggtitle(paste0("Insulin treated patients\nwith MODY probability >10%")) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ),
  
  dataset.referral_type1 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.referral_type1_with_T$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.1, "Over", "Under")) %>%
    filter(threshold == "Under") %>%
    select(M) %>%
    group_by(M) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(M = ifelse(M == 0, paste0("Non-MODY\nn=", n), paste0("MODY\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = M, values = n), colour = "white") +
    scale_fill_manual(values = c("red", "red4")) +
    theme_void() +
    ggtitle(paste0("Insulin treated patients\nwith MODY probability <10%")) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ),
  
  nrow = 1
  
) +
  patchwork::plot_annotation(
    title = "Referrals: insulin - biomarkers (10% threshold)"
  ) &
  theme(
    legend.position = "bottom"
  )


plot_referral_type1_old <- patchwork::wrap_plots(
  
  dataset.referral_type1 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.referral_type1_old$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.1, "Over", "Under")) %>%
    select(threshold) %>%
    group_by(threshold) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(threshold = ifelse(threshold == "Under", paste0("<10%\nn=", n), paste0(">10%\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = threshold, values = n), colour = "white") +
    theme_void() +
    ggtitle(paste0("Insulin treated patients\nn=", nrow(predictions_dataset.referral_type1_old))) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ),
  
  dataset.referral_type1 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.referral_type1_old$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.1, "Over", "Under")) %>%
    filter(threshold == "Over") %>%
    select(M) %>%
    group_by(M) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(M = ifelse(M == 0, paste0("Non-MODY\nn=", n), paste0("MODY\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = M, values = n), colour = "white", show.legend = TRUE) +
    scale_fill_manual(values = c("#00BFC4", "deepskyblue4")) +
    theme_void() +
    ggtitle(paste0("Insulin treated patients\nwith MODY probability >10%")) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ),
  
  dataset.referral_type1 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.referral_type1_old$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.1, "Over", "Under")) %>%
    filter(threshold == "Under") %>%
    select(M) %>%
    group_by(M) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(M = ifelse(M == 0, paste0("Non-MODY\nn=", n), paste0("MODY\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = M, values = n), colour = "white") +
    scale_fill_manual(values = c("red", "red4")) +
    theme_void() +
    ggtitle(paste0("Insulin treated patients\nwith MODY probability <10%")) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ),
  
  nrow = 1
  
) +
  patchwork::plot_annotation(
    title = "Referrals: insulin - old (10% threshold)"
  ) &
  theme(
    legend.position = "bottom"
  )



plot_referral_type2 <- patchwork::wrap_plots(
  
  dataset.referral_type2 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.referral_type2_new$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.25, "Over", "Under")) %>%
    select(threshold) %>%
    group_by(threshold) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(threshold = ifelse(threshold == "Under", paste0("<25%\nn=", n), paste0(">25%\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = threshold, values = n), colour = "white") +
    theme_void() +
    ggtitle(paste0("Non-insulin treated patients\nn=", nrow(predictions_dataset.referral_type2_new))) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ),
  
  dataset.referral_type2 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.referral_type2_new$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.25, "Over", "Under")) %>%
    filter(threshold == "Over") %>%
    select(M) %>%
    group_by(M) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(M = ifelse(M == 0, paste0("Non-MODY\nn=", n), paste0("MODY\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = M, values = n), colour = "white", show.legend = TRUE) +
    scale_fill_manual(values = c("#00BFC4", "deepskyblue4")) +
    theme_void() +
    ggtitle(paste0("Insulin treated patients\nwith MODY probability >25%")) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ),
  
  dataset.referral_type2 %>%
    select(M) %>%
    cbind(
      prop = predictions_dataset.referral_type2_new$prob
    ) %>%
    mutate(threshold = ifelse(prop > 0.25, "Over", "Under")) %>%
    filter(threshold == "Under") %>%
    select(M) %>%
    group_by(M) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    unique() %>%
    mutate(M = ifelse(M == 0, paste0("Non-MODY\nn=", n), paste0("MODY\nn=", n))) %>%
    ggplot() +
    geom_waffle(aes(fill = M, values = n), colour = "white") +
    scale_fill_manual(values = c("red", "red4")) +
    theme_void() +
    ggtitle(paste0("Non-insulin treated patients\nwith MODY probability <25%")) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(hjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ),
  
  nrow = 1
  
) +
  patchwork::plot_annotation(
    title = "Referrals: non-insulin - clinical features (25% threshold)"
  ) &
  theme(
    legend.position = "bottom"
  )




pdf("threshold_investigation.pdf", width = 12, height = 5)
plot_UNITED_type1_no_T
plot_UNITED_type1_with_T
plot_UNITED_type1_old
plot_UNITED_type2
plot_referral_type1_no_T
plot_referral_type1_with_T
plot_referral_type1_old
plot_referral_type2
dev.off()



