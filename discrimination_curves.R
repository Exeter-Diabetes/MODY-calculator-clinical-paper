#:--------------------------------------------------------
#   
# In this file we check AUROC for the MODY models
#
#:--------------------------------------------------------

# load libraries
library(tidyverse)
library(nimble)
library(pROC)

# load functions needed
source("data/create_data.R")
source("new_data_predictions/prediction_functions.R")

# load files required
predictions_dataset.UNITED_type1_no_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_no_T_full.rds")
predictions_dataset.UNITED_type1_with_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_with_T_full.rds")
# predictions_dataset.referral_type1_no_T <- readRDS("model_predictions/predictions_dataset.referral_type1_no_T_full.rds")
# predictions_dataset.referral_type1_with_T <- readRDS("model_predictions/predictions_dataset.referral_type1_with_T_full.rds")

predictions_dataset.UNITED_type2_new <- readRDS("model_predictions/predictions_dataset.UNITED_type2_new_full.rds")
# predictions_dataset.referral_type2_new <- readRDS("model_predictions/predictions_dataset.referral_type2_new_full.rds")


# load datasets

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



#:------------------------------------------------------------

# Calculate AUROC with intervals
calc_auroc <- function(data, predictions, thinning = 1) {
  
  output <- NULL
  
  for (i in seq(1, nrow(predictions), thinning)) {
    
    ## calculate auroc
    interim <- roc(data, 
                   predictions[i,],
                   plot=FALSE, print.auc=FALSE)
    
    ## append new value
    output <- c(output, as.numeric(interim$auc))
    
  }
  
  return(output)
  
}

###
## This section is thinned to help run times but also does not make a difference to full posterior values
###

## Type 1 UNITED

### No biomarker models
auc_T1D_no_T_united <- calc_auroc(dataset.UNITED_type1$M, predictions_dataset.UNITED_type1_no_T, thinning = 1000)

### Biomarker models
auc_T1D_with_T_united <- calc_auroc(dataset.UNITED_type1$M, predictions_dataset.UNITED_type1_with_T, thinning = 1000)

## Type 2 UNITED
auc_T2D_no_T_united <- calc_auroc(dataset.UNITED_type2$M, predictions_dataset.UNITED_type2_new, thinning = 1000)

# ## Type 1 referrals
# 
# ### No biomarker models
# auc_T1D_no_T_referrals <- calc_auroc(dataset.referral_type1$M, predictions_dataset.referral_type1_no_T, thinning = 1000)
# 
# ### Biomarker models
# auc_T1D_with_T_referrals <- calc_auroc(dataset.referral_type1$M, predictions_dataset.referral_type1_with_T, thinning = 1000)
# 
# ## Type 2 referrals
# auc_T2D_no_T_referrals <- calc_auroc(dataset.referral_type2$M, predictions_dataset.referral_type2_new, thinning = 1000)


#:------------------------------------------------------------

# Calculate ROC with intervals
calc_roc <- function(data, predictions, thinning = 1) {
  
  output <- NULL
  
  for (i in seq(1, nrow(predictions), thinning)) {
    
    ## calculate ROC
    interim <- pROC::roc(response = data, predictor = predictions[i,]) %>%
      magrittr::extract(c(2:3)) %>%
      as.data.frame() %>%
      mutate(iteration = paste0(i))
    
    output <- rbind(output, interim)
    
  }
  
  return(output)
  
}


## Type 1 UNITED

### No biomarker models
roc_T1D_no_T_united <- calc_roc(dataset.UNITED_type1$M, predictions_dataset.UNITED_type1_no_T, thinning = 100)

# plot for ROC with grey being iterations, black being the ROC for average prediction
plot_roc_T1D_no_T_united <- ggplot() +
  ## all iterations
  geom_path(
    data = roc_T1D_no_T_united,
    aes(x = 1-specificities, y= sensitivities, group = iteration), colour = "grey"
  ) +
  ## average predictions
  geom_path(
    data = pROC::roc(response = dataset.UNITED_type1$M, predictor = colMeans(predictions_dataset.UNITED_type1_no_T)) %>%
      magrittr::extract(c(2:3)) %>%
      as.data.frame(),
    aes(x = 1-specificities, y= sensitivities), colour = "black"
  )

### Biomarker models
roc_T1D_with_T_united <- calc_roc(dataset.UNITED_type1$M, predictions_dataset.UNITED_type1_with_T, thinning = 100)

# plot for ROC with grey being iterations, black being the ROC for average prediction
plot_roc_T1D_with_T_united <- ggplot() +
  ## all iterations
  geom_path(
    data = roc_T1D_with_T_united,
    aes(x = 1-specificities, y= sensitivities, group = iteration), colour = "grey"
  ) +
  ## average predictions
  geom_path(
    data = pROC::roc(response = dataset.UNITED_type1$M, predictor = colMeans(predictions_dataset.UNITED_type1_with_T)) %>%
      magrittr::extract(c(2:3)) %>%
      as.data.frame(),
    aes(x = 1-specificities, y= sensitivities), colour = "black"
  )

## Type 2 UNITED
roc_T2D_new_united <- calc_roc(dataset.UNITED_type2$M, predictions_dataset.UNITED_type2_new, thinning = 100)

# plot for ROC with grey being iterations, black being the ROC for average prediction
plot_roc_T2D_new_united <- ggplot() +
  ## all iterations
  geom_path(
    data = roc_T2D_new_united,
    aes(x = 1-specificities, y= sensitivities, group = iteration), colour = "grey"
  ) +
  ## average predictions
  geom_path(
    data = pROC::roc(response = dataset.UNITED_type2$M, predictor = colMeans(predictions_dataset.UNITED_type2_new)) %>%
      magrittr::extract(c(2:3)) %>%
      as.data.frame(),
    aes(x = 1-specificities, y= sensitivities), colour = "black"
  )


# ## Type 1 referrals
# 
# ### No biomarker models
# roc_T1D_no_T_referral <- calc_roc(dataset.referral_type1$M, predictions_dataset.referral_type1_no_T, thinning = 1000)
# 
# # plot for ROC with grey being iterations, black being the ROC for average prediction
# plot_roc_T1D_no_T_referral <- ggplot() +
#   ## all iterations
#   geom_path(
#     data = roc_T1D_no_T_referral,
#     aes(x = 1-specificities, y= sensitivities, group = iteration), colour = "grey"
#   ) +
#   ## average predictions
#   geom_path(
#     data = pROC::roc(response = dataset.referral_type1$M, predictor = colMeans(predictions_dataset.referral_type1_no_T)) %>%
#       magrittr::extract(c(2:3)) %>%
#       as.data.frame(),
#     aes(x = 1-specificities, y= sensitivities), colour = "black"
#   )
# 
# ### Biomarker models
# roc_T1D_with_T_referral <- calc_roc(dataset.referral_type1$M, predictions_dataset.referral_type1_with_T, thinning = 1000)
# 
# # plot for ROC with grey being iterations, black being the ROC for average prediction
# plot_roc_T1D_with_T_referral <- ggplot() +
#   ## all iterations
#   geom_path(
#     data = roc_T1D_with_T_referral,
#     aes(x = 1-specificities, y= sensitivities, group = iteration), colour = "grey"
#   ) +
#   ## average predictions
#   geom_path(
#     data = pROC::roc(response = dataset.referral_type1$M, predictor = colMeans(predictions_dataset.referral_type1_with_T)) %>%
#       magrittr::extract(c(2:3)) %>%
#       as.data.frame(),
#     aes(x = 1-specificities, y= sensitivities), colour = "black"
#   )
# 
# ## Type 2 referral
# roc_T2D_new_referral <- calc_roc(dataset.referral_type2$M, predictions_dataset.referral_type2_new, thinning = 1000)
# 
# # plot for ROC with grey being iterations, black being the ROC for average prediction
# plot_roc_T2D_new_referral <- ggplot() +
#   ## all iterations
#   geom_path(
#     data = roc_T2D_new_referral,
#     aes(x = 1-specificities, y= sensitivities, group = iteration), colour = "grey"
#   ) +
#   ## average predictions
#   geom_path(
#     data = pROC::roc(response = dataset.referral_type2$M, predictor = colMeans(predictions_dataset.referral_type2_new)) %>%
#       magrittr::extract(c(2:3)) %>%
#       as.data.frame(),
#     aes(x = 1-specificities, y= sensitivities), colour = "black"
#   )




#:------------------------------------------------------------


# Calculate Precision-recall with intervals
calc_prec_recal_curve <- function(data, predictions, thinning = 1) {
  
  output <- NULL
  
  for (i in seq(1, nrow(predictions), thinning)) {
    
    ## calculate ROC
    interim <- pROC::coords(pROC::roc(response = data, predictor = predictions[i,]), ret = c("precision", "recall")) %>%
      as.data.frame() %>%
      mutate(iteration = paste0(i))
    
    output <- rbind(output, interim)
    
  }
  
  return(output)
  
}

## Type 1 UNITED

### No biomarker models
prec_recal_T1D_no_T_united <- calc_prec_recal_curve(dataset.UNITED_type1$M, predictions_dataset.UNITED_type1_no_T, thinning = 1000)

# plot for ROC with grey being iterations, black being the ROC for average prediction
plot_prec_recal_T1D_no_T_united <- ggplot() +
  ## all iterations
  geom_path(
    data = prec_recal_T1D_no_T_united,
    aes(x = recall, y= precision, group = iteration), colour = "grey"
  ) +
  ## average predictions
  geom_path(
    data = pROC::coords(pROC::roc(response = dataset.UNITED_type1$M, predictor = colMeans(predictions_dataset.UNITED_type1_no_T)), ret = c("precision", "recall")) %>%
      as.data.frame(),
    aes(x = recall, y= precision), colour = "black"
  )

### Biomarker models
prec_recal_T1D_with_T_united <- calc_prec_recal_curve(dataset.UNITED_type1$M, predictions_dataset.UNITED_type1_with_T, thinning = 1000)

# plot for ROC with grey being iterations, black being the ROC for average prediction
plot_prec_recal_T1D_with_T_united <- ggplot() +
  ## all iterations
  geom_path(
    data = prec_recal_T1D_with_T_united,
    aes(x = recall, y= precision, group = iteration), colour = "grey"
  ) +
  ## average predictions
  geom_path(
    data = pROC::coords(pROC::roc(response = dataset.UNITED_type1$M, predictor = colMeans(predictions_dataset.UNITED_type1_with_T)), ret = c("precision", "recall")) %>%
      as.data.frame(),
    aes(x = recall, y= precision), colour = "black"
  )

## Type 2 UNITED
prec_recal_T2D_new_united <- calc_prec_recal_curve(dataset.UNITED_type2$M, predictions_dataset.UNITED_type2_new, thinning = 1000)

# plot for ROC with grey being iterations, black being the ROC for average prediction
plot_prec_recal_T2D_new_united <- ggplot() +
  ## all iterations
  geom_path(
    data = prec_recal_T2D_new_united,
    aes(x = recall, y= precision, group = iteration), colour = "grey"
  ) +
  ## average predictions
  geom_path(
    data = pROC::coords(pROC::roc(response = dataset.UNITED_type2$M, predictor = colMeans(predictions_dataset.UNITED_type2_new)), ret = c("precision", "recall")) %>%
      as.data.frame(),
    aes(x = recall, y= precision), colour = "black"
  )


# ## Type 1 referral
# 
# ### No biomarker models
# prec_recal_T1D_no_T_referral <- calc_prec_recal_curve(dataset.referral_type1$M, predictions_dataset.referral_type1_no_T, thinning = 1000)
# 
# # plot for ROC with grey being iterations, black being the ROC for average prediction
# plot_prec_recal_T1D_no_T_referral <- ggplot() +
#   ## all iterations
#   geom_path(
#     data = prec_recal_T1D_no_T_referral,
#     aes(x = recall, y= precision, group = iteration), colour = "grey"
#   ) +
#   ## average predictions
#   geom_path(
#     data = pROC::coords(pROC::roc(response = dataset.referral_type1$M, predictor = colMeans(predictions_dataset.referral_type1_no_T)), ret = c("precision", "recall")) %>%
#       as.data.frame(),
#     aes(x = recall, y= precision), colour = "black"
#   )
# 
# ### Biomarker models
# prec_recal_T1D_with_T_referral <- calc_prec_recal_curve(dataset.referral_type1$M, predictions_dataset.referral_type1_with_T, thinning = 1000)
# 
# # plot for ROC with grey being iterations, black being the ROC for average prediction
# plot_prec_recal_T1D_with_T_referral <- ggplot() +
#   ## all iterations
#   geom_path(
#     data = prec_recal_T1D_with_T_referral,
#     aes(x = recall, y= precision, group = iteration), colour = "grey"
#   ) +
#   ## average predictions
#   geom_path(
#     data = pROC::coords(pROC::roc(response = dataset.referral_type1$M, predictor = colMeans(predictions_dataset.referral_type1_with_T)), ret = c("precision", "recall")) %>%
#       as.data.frame(),
#     aes(x = recall, y= precision), colour = "black"
#   )
# 
# ## Type 2 referral
# prec_recal_T2D_new_referral <- calc_prec_recal_curve(dataset.referral_type2$M, predictions_dataset.referral_type2_new, thinning = 1000)
# 
# # plot for ROC with grey being iterations, black being the ROC for average prediction
# plot_prec_recal_T2D_new_referral <- ggplot() +
#   ## all iterations
#   geom_path(
#     data = prec_recal_T2D_new_referral,
#     aes(x = recall, y= precision, group = iteration), colour = "grey"
#   ) +
#   ## average predictions
#   geom_path(
#     data = pROC::coords(pROC::roc(response = dataset.referral_type2$M, predictor = colMeans(predictions_dataset.referral_type2_new)), ret = c("precision", "recall")) %>%
#       as.data.frame(),
#     aes(x = recall, y= precision), colour = "black"
#   )




#:--------------------------------------------------

# Boxplot and roc curves

roc_curves <- data.frame(prob = colMeans(predictions_dataset.UNITED_type1_with_T)) %>%
  cbind(Mody = dataset.UNITED_type1$M) %>%
  pROC::roc(response = Mody, predictor = prob) %>%
  magrittr::extract(2:3) %>%
  as.data.frame() %>%
  mutate(
    ROCAUC =  unname(data.frame(prob = colMeans(predictions_dataset.UNITED_type1_with_T)) %>%
                       cbind(Mody = dataset.UNITED_type1$M) %>%
                       pROC::roc(response = Mody, predictor = prob) %>%
                       magrittr::extract(c(9)) %>%
                       unlist()),
    mean = mean(colMeans(predictions_dataset.UNITED_type1_with_T), na.rm = TRUE)
  ) %>%
  mutate(Dataset = "UNITED", Model = "Type 1", Calculator = "Biomarkers") %>%
  rbind(
    data.frame(prob = colMeans(predictions_dataset.UNITED_type1_no_T)) %>%
      cbind(Mody = dataset.UNITED_type1$M) %>%
      pROC::roc(response = Mody, predictor = prob) %>%
      magrittr::extract(2:3) %>%
      as.data.frame() %>%
      mutate(
        ROCAUC = unname(data.frame(prob = colMeans(predictions_dataset.UNITED_type1_no_T)) %>%
                          cbind(Mody = dataset.UNITED_type1$M) %>%
                          pROC::roc(response = Mody, predictor = prob) %>%
                          magrittr::extract(c(9)) %>%
                          unlist()),
        mean = mean(colMeans(predictions_dataset.UNITED_type1_no_T), na.rm = TRUE)
      ) %>%
      mutate(Dataset = "UNITED", Model = "Type 1", Calculator = "No Biomarkers"), 
    data.frame(prob = colMeans(predictions_dataset.UNITED_type2_new)) %>%
      cbind(Mody = dataset.UNITED_type2$M) %>%
      pROC::roc(response = Mody, predictor = prob) %>%
      magrittr::extract(2:3) %>%
      as.data.frame() %>%
      mutate(
        ROCAUC = unname(data.frame(prob = colMeans(predictions_dataset.UNITED_type2_new)) %>%
                          cbind(Mody = dataset.UNITED_type2$M) %>%
                          pROC::roc(response = Mody, predictor = prob) %>%
                          magrittr::extract(c(9)) %>%
                          unlist()),
        mean = mean(colMeans(predictions_dataset.UNITED_type2_new), na.rm = TRUE)
      ) %>%
      mutate(Dataset = "UNITED", Model = "Type 2", Calculator = " ")
    # , data.frame(prob = colMeans(predictions_dataset.referral_type1_with_T)) %>%
    #   cbind(Mody = dataset.referral_type1$M) %>%
    #   pROC::roc(response = Mody, predictor = prob) %>%
    #   magrittr::extract(2:3) %>%
    #   as.data.frame() %>%
    #   mutate(
    #     ROCAUC = unname(data.frame(prob = colMeans(predictions_dataset.referral_type1_with_T)) %>%
    #       cbind(Mody = dataset.referral_type1$M) %>%
    #       pROC::roc(response = Mody, predictor = prob) %>%
    #         magrittr::extract(c(9)) %>%
    #         unlist()),
    #     mean = mean(colMeans(predictions_dataset.referral_type1_with_T), na.rm = TRUE)
    #   ) %>%
    #   mutate(Dataset = "Referral", Model = "Type 1", Calculator = "Biomarkers"),
    # data.frame(prob = colMeans(predictions_dataset.referral_type1_no_T)) %>%
    #   cbind(Mody = dataset.referral_type1$M) %>%
    #   pROC::roc(response = Mody, predictor = prob) %>%
    #   magrittr::extract(2:3) %>%
    #   as.data.frame() %>%
    #   mutate(
    #     ROCAUC = unname(data.frame(prob = colMeans(predictions_dataset.referral_type1_no_T)) %>%
    #       cbind(Mody = dataset.referral_type1$M) %>%
    #       pROC::roc(response = Mody, predictor = prob) %>%
    #         magrittr::extract(c(9)) %>%
    #         unlist()),
    #     mean = mean(colMeans(predictions_dataset.referral_type1_no_T), na.rm = TRUE)
    #   ) %>%
    #   mutate(Dataset = "Referral", Model = "Type 1", Calculator = "No Biomarkers"),
    # data.frame(prob = colMeans(predictions_dataset.referral_type2_new)) %>%
    #   cbind(Mody = dataset.referral_type2$M) %>%
    #   pROC::roc(response = Mody, predictor = prob) %>%
    #   magrittr::extract(2:3) %>%
    #   as.data.frame() %>%
    #   mutate(
    #     ROCAUC = unname(data.frame(prob = colMeans(predictions_dataset.referral_type2_new)) %>%
    #       cbind(Mody = dataset.referral_type2$M) %>%
    #       pROC::roc(response = Mody, predictor = prob) %>%
    #         magrittr::extract(c(9)) %>%
    #         unlist()),
    #     mean = mean(colMeans(predictions_dataset.referral_type2_new), na.rm = TRUE)
    #   ) %>%
    #   mutate(Dataset = "Referral", Model = "Type 2", Calculator = " ")
  ) %>%
  rename("auc" = "ROCAUC")


dat_text <- roc_curves %>%
  select(-sensitivities, -specificities) %>%
  distinct()

dat_text$ROCAUC <- unlist(dat_text$ROCAUC)

dat_text <- dat_text %>%
  mutate(
    ROCAUC = paste0(" AUC:", signif(ROCAUC, 2), " "),
    mean = paste0("Mean prob:", signif(mean, 2)*100, "%")
  )



plot_prob_boxplot_rocs_united <- patchwork::wrap_plots(
  
  # Panel A - UNITED insulin-treated
  
  patchwork::wrap_plots(
    
    # Boxplots
    dataset.UNITED_type1 %>%
      select(M) %>%
      rename("Mody" = "M") %>%
      cbind(
        prob_with = colMeans(predictions_dataset.UNITED_type1_with_T),
        prob_without = colMeans(predictions_dataset.UNITED_type1_no_T)
      ) %>%
      gather("key", "Probability", -Mody) %>% 
      mutate(
        Mody = factor(Mody, levels = c(0, 1), labels = c("Negative", "Positive")),
        key = factor(key, levels = c("prob_without", "prob_with"), labels = c("No Biomarkers", "Biomarkers"))
      ) %>%
      ggplot() +
      geom_boxplot(aes(y = Probability, x = Mody), colour = c("black", "white", "black", "white"), alpha = c(1, 0, 1, 0)) +
      geom_point(aes(y = Probability, x = Mody, colour = Mody, alpha = Mody)) +
      scale_y_continuous(labels = scales::percent) +
      scale_colour_manual(values = c("white", "black")) +
      scale_alpha_manual(values = c(0, 1)) +
      facet_wrap(~key) +
      theme_bw() +
      theme(
        legend.position = "none"
      ),
    roc_curves %>%
      filter(Dataset == "UNITED" & Model == "Type 1") %>%
      mutate(
        Calculator = factor(Calculator, levels = c("Biomarkers", "No Biomarkers")),
        iteration = 0
      ) %>%
      ggplot(aes(x = 1- specificities, y = sensitivities)) +
      geom_path(
        data = roc_T1D_no_T_united %>%
          mutate(
            Calculator = "No Biomarkers"
          ) %>%
          rbind(
            roc_T1D_with_T_united %>%
              mutate(Calculator = "Biomarkers")
          ),
        aes(group = iteration), colour = "grey"
      ) + 
      geom_path() +
      theme_bw() +
      facet_grid(~factor(Calculator, levels = c("No Biomarkers", "Biomarkers")), scales = "free",) +
      scale_y_continuous("Sensitivity", labels = scales::percent) +
      scale_x_continuous("1- Specificity", labels = scales::percent) +
      theme_bw() +
      geom_label(
        data = dat_text %>%
          filter(Dataset == "UNITED" & Model == "Type 1"),
        mapping = aes(x = -Inf, y = -Inf, label = ROCAUC),
        size = 7,
        label.size = NA,
        hjust = -0.3,
        vjust = -0.5
      ) +
      theme(
        panel.spacing.x = unit(1.5, "lines")
      ), 
    
    ncol = 2, nrow = 1
  ),
  
  
  # Panel B - UNITED non-insulin-treated
  
  patchwork::wrap_plots(
    
    # Boxplots
    dataset.UNITED_type2 %>%
      select(M) %>%
      rename("Mody" = "M") %>%
      cbind(
        Probability = colMeans(predictions_dataset.UNITED_type2_new),
        key = ""
      ) %>%
      mutate(
        Mody = factor(Mody, levels = c(0, 1), labels = c("Negative", "Positive"))
      ) %>%
      ggplot() +
      geom_boxplot(aes(y = Probability, x = Mody)) +
      scale_y_continuous(labels = scales::percent) +
      theme_bw(),
    
    roc_curves %>%
      filter(Dataset == "UNITED" & Model == "Type 2") %>%
      ggplot(aes(x = 1- specificities, y = sensitivities)) +
      geom_path(
        data = roc_T2D_new_united %>%
          mutate(
            Calculator = "No Biomarkers"
          ),
        aes(group = iteration), colour = "grey"
      ) +
      geom_path() +
      theme_bw() +
      scale_y_continuous("Sensitivity", labels = scales::percent) +
      scale_x_continuous("1- Specificity", labels = scales::percent) +
      theme_bw() +
      geom_label(
        data = dat_text %>%
          filter(Dataset == "UNITED" & Model == "Type 2"),
        mapping = aes(x = -Inf, y = -Inf, label = ROCAUC),
        size = 7,
        label.size = NA,
        hjust = -0.4,
        vjust = -0.5
      ),
    
    
    ncol = 2, nrow = 1
    
  ),
  
  
  ncol = 1
  
  
) + patchwork::plot_annotation(tag_levels = list(c("A.1", "A.2", "B.1", "B.2"))) &
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 15)
  )


#:-------------------------------------------------------------
# Making plots
pdf("figures/united_boxplot_roc_thin_100.pdf", width = 11, height = 9)
plot_prob_boxplot_rocs_united
dev.off()





plot_prob_boxplot_rocs_referral <- patchwork::wrap_plots(
  
  # Panel A - Referral insulin-treated
  
  patchwork::wrap_plots(
    
    # Boxplots
    dataset.referral_type1 %>%
      select(M) %>%
      rename("Mody" = "M") %>%
      cbind(
        prob_with = colMeans(predictions_dataset.referral_type1_with_T),
        prob_without = colMeans(predictions_dataset.referral_type1_no_T)
      ) %>%
      gather("key", "Probability", -Mody) %>% 
      mutate(
        Mody = factor(Mody, levels = c(0, 1), labels = c("Negative", "Positive")),
        key = factor(key, levels = c("prob_without", "prob_with"), labels = c("No Biomarkers", "Biomarkers"))
      ) %>%
      ggplot() +
      geom_boxplot(aes(y = Probability, x = Mody)) +
      scale_y_continuous(labels = scales::percent) +
      facet_wrap(~key) +
      theme_bw() +
      theme(
        legend.position = "none"
      ),
    roc_curves %>%
      filter(Dataset == "Referral" & Model == "Type 1") %>%
      mutate(
        Calculator = factor(Calculator, levels = c("No Biomarkers", "Biomarkers"))
      ) %>%
      ggplot(aes(x = 1- specificities, y = sensitivities)) +
      geom_path() +
      theme_bw() +
      facet_grid(~Calculator, scales = "free") +
      scale_y_continuous("Sensitivity", labels = scales::percent) +
      scale_x_continuous("1- Specificity", labels = scales::percent) +
      theme_bw() +
      geom_label(
        data = dat_text %>%
          filter(Dataset == "Referral" & Model == "Type 1"),
        mapping = aes(x = -Inf, y = -Inf, label = ROCAUC),
        size = 7,
        label.size = NA,
        hjust = -0.4,
        vjust = -0.5
      ), 
    
    ncol = 2, nrow = 1
  ),
  
  
  # Panel B - Referral non-insulin-treated
  
  patchwork::wrap_plots(
    
    # Boxplots
    dataset.referral_type2 %>%
      select(M) %>%
      rename("Mody" = "M") %>%
      cbind(
        Probability = colMeans(predictions_dataset.referral_type2_new),
        key = ""
      ) %>%
      mutate(
        Mody = factor(Mody, levels = c(0, 1), labels = c("Negative", "Positive"))
      ) %>%
      ggplot() +
      geom_boxplot(aes(y = Probability, x = Mody)) +
      scale_y_continuous(labels = scales::percent) +
      facet_wrap(~key) +
      theme_bw(),
    
    roc_curves %>%
      filter(Dataset == "Referral" & Model == "Type 2") %>%
      ggplot(aes(x = 1- specificities, y = sensitivities)) +
      geom_path() +
      theme_bw() +
      facet_grid(~Calculator, scales = "free") +
      scale_y_continuous("Sensitivity", labels = scales::percent) +
      scale_x_continuous("1- Specificity", labels = scales::percent) +
      theme_bw() +
      geom_label(
        data = dat_text %>%
          filter(Dataset == "Referral" & Model == "Type 2"),
        mapping = aes(x = -Inf, y = -Inf, label = ROCAUC),
        size = 7,
        label.size = NA,
        hjust = -0.4,
        vjust = -0.5
      ),
    
    
    ncol = 2, nrow = 1
    
  ),
  
  
  ncol = 1
  
  
) + patchwork::plot_annotation(tag_levels = list(c("A.1", "A.2", "B.1", "B.2")))








