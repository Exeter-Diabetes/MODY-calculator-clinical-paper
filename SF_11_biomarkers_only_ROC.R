#################################################################################

#Supplementary Figure 11

#ROC curve of using biomarker status to predict MODY probability 
#for early-insulin-treated patients

######################################################################################
setwd("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper")
#load libraries --------------------------------------------------------------------------------
library(tidyverse)
library(pROC)
library(PRROC)
library(writexl)

#load functions
source("data/create_data.R")
## Load data  --------------------------------------------------------------------------
dataset.UNITED_type1_all_genes <- create_data(dataset = "united t1d", commonmody = FALSE, id = TRUE) %>%
  ## if MODY testing missing, change to 0
  mutate(M = ifelse(is.na(M), 0, M))

#load Early-insulin-treated predictions
predictions_dataset.UNITED_type1_all_genes_with_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_all_genes_with_T.rds")
UNITED_type1 <- cbind(dataset.UNITED_type1_all_genes, predictions_dataset.UNITED_type1_all_genes_with_T)

roc(UNITED_type1$M, UNITED_type1$T, plot = TRUE, print.thres = "best", print.auc = TRUE, ci = TRUE)

# Calculate AUROC with intervals
calc_auroc <- function(data, predictions, thinning = 100) {
  
  # output file
  output <- NULL
  
  # sequence to iterate through
  sequence_list <- seq(1, nrow(predictions), thinning)
  
  for (i in 1:length(sequence_list)) {
    
    # print the current iteration
    if (i %% 1000 == 0) {
      print(paste(i, "out of", length(sequence_list)))
    }
    
    interim <- roc(data, 
                   predictions[sequence_list[i],],
                   levels=c(0,1), direction = "<", plot=FALSE, print.auc=FALSE)
    
    ## append new value
    output <- c(output, as.numeric(interim$auc))
    
  }
  
  return(output)
  
}

###
## This section is thinned to help run times but also does not make much difference to full posterior values
###

## Type 1 UNITED

### Biomarker models
auc_T1D_with_T_united_all_genes <- calc_auroc(UNITED_type1$M, predictions_dataset.UNITED_type1_all_genes_with_T_full, thinning = 10)
# quantile(auc_T1D_with_T_united_all_genes, probs = c(0.025, 0.5, 0.975)) # thinning = 10
# 2.5%       50%     97.5% 
# 0.9487113 0.9768041 0.9779210 

roc_curves <- data.frame(T = UNITED_type1$T) %>%
  cbind(Mody = UNITED_type1$M) %>%
  pROC::roc(response = Mody, predictor = T) %>%
  magrittr::extract(2:3) %>%
  as.data.frame() %>%
  mutate(
    auc =  unname(data.frame(prob = UNITED_type1$T) %>%
                    cbind(Mody = UNITED_type1$M) %>%
                    pROC::roc(response = Mody, predictor = T) %>%
                    magrittr::extract(c(9)) %>%
                    unlist()),
    auc_low = quantile(auc_T1D_with_T_united_all_genes, probs = c(0.025)),
    auc_high = quantile(auc_T1D_with_T_united_all_genes, probs= c(0.975)),
    mean = mean(UNITED_type1$T, na.rm = TRUE)
  ) %>%
  mutate(Dataset = "UNITED", Model = "Type 1", Calculator = "Biomarkers") 


dat_text <- roc_curves %>%
  select(-sensitivities, -specificities) %>%
  distinct() %>%
  mutate(
    auc_full = paste0("AUC: ", signif(auc, 2), " [", signif(auc_low, 2), "-", signif(auc_high, 2), "]"),
    mean = paste0("Mean prob:", signif(mean, 2)*100, "%"),
    Calculator = factor(Calculator, levels = c("Biomarkers"), labels = c("Biomarkers"))
  )

roc_curves %>%
  filter(Dataset == "UNITED" & Model == "Type 1" & Calculator == "Biomarkers") %>%
  mutate(
    Calculator = factor(Calculator, levels = c("Biomarkers"), labels = c("Biomarkers")),
    iteration = 0
  ) %>%
  ggplot(aes(x = 1- specificities, y = sensitivities)) +
  geom_path() +
  theme_bw() +
  # facet_grid(~factor(Calculator, levels = c("Clinical features", "Clinical features and biomarkers"), labels = c("Early-insulin-treated: clinical features, no biomarkers", "Early-insulin-treated: clinical features, with biomarkers")), scales = "free",labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
  scale_y_continuous("Sensitivity", labels = scales::percent) +
  scale_x_continuous("1- Specificity", labels = scales::percent) +
  theme_bw() +
  geom_label(
    data = dat_text %>%
      filter(Dataset == "UNITED" & Model == "Type 1" & Calculator == "Biomarkers"),
    mapping = aes(x = 0.55, y = 0.1, label = auc_full, hjust = "center"),
    size = 7,
    label.r = unit(0, "pt"),
    label.padding=unit(0.4, "lines")
  )


t.test(auc_T1D_with_T_united_all_genes, auc_T1D_no_T_united_all_genes)
