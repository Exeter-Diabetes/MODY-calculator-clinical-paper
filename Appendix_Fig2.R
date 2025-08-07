#######################################################################################
#Appendix Figure 

#Shows precision-recall curves for UNITED internal validation
###########################################################################################
# load libraries
library(tidyverse)
library(nimble)
library(pROC)
library(PRROC)


# load functions needed
source("Data/create_data.R")
source("New_Data_Predictions/prediction_functions.R")

# load datasets
## Load population representative dataset
dataset_UNITED_type1 <- create_data(
  dataset = "united t1d", 
  commonmody = FALSE, 
  id = TRUE) %>%
  
  ## if MODY testing missing, change to 0
  mutate(M = ifelse(is.na(M), 0, M))

dataset_UNITED_type2 <- create_data(
  dataset = "united t2d", 
  commonmody = FALSE, 
  id = TRUE)

# load files required
predictions_UNITED_type1_no_T_full <- readRDS(
  "Model_Predictions/predictions_dataset.UNITED_type1_all_genes_no_T_full.rds")
predictions_UNITED_type1_no_T_full <- predictions_UNITED_type1_no_T_full[, as.character(c(dataset_UNITED_type1$id))]
predictions_UNITED_type1_with_T_full <- readRDS(
  "Model_Predictions/predictions_dataset.UNITED_type1_all_genes_with_T_full.rds")
predictions_UNITED_type1_with_T_full <- predictions_UNITED_type1_with_T_full[, as.character(c(dataset_UNITED_type1$id))]
predictions_UNITED_type2_new_full <- readRDS(
  "Model_Predictions/predictions_dataset.UNITED_type2_all_genes_new_full.rds")
predictions_UNITED_type2_new_full <- predictions_UNITED_type2_new_full[, as.character(c(dataset_UNITED_type2$id))]

predictions_UNITED_type1_no_T <- readRDS(
  "Model_Predictions/predictions_dataset.UNITED_type1_all_genes_no_T.rds") %>% 
  as.data.frame() %>%
  { rownames(.) <- NULL; . } %>%
  column_to_rownames(var = "id")
predictions_UNITED_type1_no_T <- predictions_UNITED_type1_no_T[as.character(c(dataset_UNITED_type1$id)), ]
predictions_UNITED_type1_with_T <- readRDS(
  "Model_Predictions/predictions_dataset.UNITED_type1_all_genes_with_T.rds") %>% 
  as.data.frame() %>%
  { rownames(.) <- NULL; . } %>%
  column_to_rownames(var = "id")
predictions_UNITED_type1_with_T <- predictions_UNITED_type1_with_T[as.character(c(dataset_UNITED_type1$id)), ]
predictions_UNITED_type2_new <- readRDS(
  "Model_Predictions/predictions_dataset.UNITED_type2_all_genes_new.rds") %>% 
  as.data.frame() %>%
  { rownames(.) <- NULL; . } %>%
  column_to_rownames(var = "id")
predictions_UNITED_type2_new <- predictions_UNITED_type2_new[as.character(c(dataset_UNITED_type2 <- create_data(dataset = "united t2d", commonmody = FALSE, id = TRUE)
                                                                            $id)), ]


#:------------------------------------------------------------
# Calculate AUROC with intervals
calc_auc_pr <- function(data, predictions, class1 = 1, thinning = 100) {
  
  # output file
  output <- NULL
  
  # sequence to iterate through
  sequence_list <- seq(1, nrow(predictions), thinning)
  
  for (i in 1:length(sequence_list)) {
    
    # print the current iteration
    if (i %% 1000 == 0) {
      print(paste(i, "out of", length(sequence_list)))
    }
    
    interim <- data.frame(
      resp = data,
      pred = predictions[sequence_list[i],]
    )
    
    ## calculate auc pr
    interim_curve <- pr.curve(
      scores.class1 = interim %>% filter(resp != class1) %>% select(pred) %>% unlist(), 
      scores.class0 = interim %>% filter(resp == class1) %>% select(pred) %>% unlist(), 
      curve = FALSE)
    
    ## append new value
    output <- c(output, as.numeric(interim_curve$auc.integral))
    
  }
  
  return(output)
  
}

## Type 1 UNITED

### No biomarker models
pr_auc_T1D_no_T_united <- calc_auc_pr(
  dataset_UNITED_type1$M, 
  predictions_UNITED_type1_no_T_full,
  thinning = 10)
# quantile(pr_auc_T1D_no_T_united_all_genes, probs = c(0.025, 0.5, 0.975)) # thinning = 10
# 2.5%        50%      97.5% 
# 0.01548813 0.02644908 0.03477454 

### Biomarker models
pr_auc_T1D_with_T_united <- calc_auc_pr(
  dataset_UNITED_type1$M, 
  predictions_UNITED_type1_with_T_full,
  thinning = 10)
# quantile(pr_auc_T1D_with_T_united_all_genes, probs = c(0.025, 0.5, 0.975)) # thinning = 10
# 2.5%        50%      97.5% 
# 0.09617254 0.19216510 0.22411871 

## Type 2 UNITED
pr_auc_T2D_new_united_all <- calc_auc_pr(
  dataset_UNITED_type2$M, 
  predictions_UNITED_type2_new_full,
  thinning = 10)
# quantile(pr_auc_T2D_new_united_all_genes, probs = c(0.025, 0.5, 0.975)) # thinning = 10
# 2.5%       50%     97.5% 
# 0.5980896 0.6291710 0.6611516 

prec_recal_curves <- data.frame(
  prob = colMeans(predictions_UNITED_type1_with_T_full)) %>%
  cbind(Mody = dataset_UNITED_type1$M) %>%
  pROC::roc(response = Mody, predictor = prob) %>%
  pROC::coords(ret = c("precision", "recall")) %>%
  as.data.frame() %>%
  mutate(
    auc = calc_auc_pr(
      dataset_UNITED_type1$M, 
      t(as.data.frame(
        colMeans(predictions_UNITED_type1_with_T_full))), 
      thinning = 1),
    auc_low = quantile(pr_auc_T1D_with_T_united, probs = c(0.025)),
    auc_high = quantile(pr_auc_T1D_with_T_united, probs= c(0.975)),
    mean = mean(colMeans(predictions_UNITED_type1_with_T_full), 
                na.rm = TRUE)
  ) %>%
  mutate(Dataset = "UNITED", Model = "Type 1", Calculator = "Biomarkers") %>%
  rbind(
    data.frame(prob = colMeans(predictions_UNITED_type1_no_T_full)) %>%
      cbind(Mody = dataset_UNITED_type1$M) %>%
      pROC::roc(response = Mody, predictor = prob) %>%
      pROC::coords(ret = c("precision", "recall")) %>%
      as.data.frame() %>%
      mutate(
        auc = calc_auc_pr(
          dataset_UNITED_type1$M, 
          t(as.data.frame(colMeans(
            predictions_UNITED_type1_no_T_full))), 
          thinning = 1),
        auc_low = quantile(pr_auc_T1D_no_T_united, probs = c(0.025)),
        auc_high = quantile(pr_auc_T1D_no_T_united, probs= c(0.975)),
        mean = mean(colMeans(
          predictions_UNITED_type1_no_T_full), 
          na.rm = TRUE)
      ) %>%
      mutate(Dataset = "UNITED", Model = "Type 1", Calculator = "No Biomarkers"), 
    data.frame(prob = colMeans(predictions_UNITED_type2_new_full)) %>%
      cbind(Mody = dataset_UNITED_type2$M) %>%
      pROC::roc(response = Mody, predictor = prob) %>%
      pROC::coords(ret = c("precision", "recall")) %>%
      as.data.frame() %>%
      mutate(
        auc = calc_auc_pr(
          dataset_UNITED_type2$M, 
          t(as.data.frame(colMeans(
            predictions_UNITED_type2_new_full))), 
          thinning = 1),
        auc_low = quantile(pr_auc_T2D_new_united, probs = c(0.025)),
        auc_high = quantile(pr_auc_T2D_new_united, probs= c(0.975)),
        mean = mean(colMeans(predictions_UNITED_type2_new_full), na.rm = TRUE)
      ) %>%
      mutate(Dataset = "UNITED", Model = "Type 2", Calculator = "No Biomarkers")
  )


dat_text <- prec_recal_curves %>%
  select(-precision, -recall) %>%
  distinct() %>%
  mutate(
    auc_full = paste0("AUC: ", signif(auc, 2), 
                      " [", signif(auc_low, 2), 
                      "-", signif(auc_high, 2), "]"),
    mean = paste0("Mean prob:", signif(mean, 2)*100, "%"),
    Calculator = factor(
      Calculator, 
      levels = c("Biomarkers", "No Biomarkers"), 
      labels = c("Clinical features and biomarkers", "Clinical features"))
  )



plot_prob_prec_recal_united <- patchwork::wrap_plots(
  
  # ROC T1D
  patchwork::free(
    
    prec_recal_curves %>%
      filter(Dataset == "UNITED" & Model == "Type 1") %>%
      mutate(
        Calculator = factor(
          Calculator, 
          levels = c("No Biomarkers", "Biomarkers"), 
          labels = c("Clinical features", "Clinical features and biomarkers")),
        iteration = 0
      ) %>%
      ggplot(aes(x = recall, y = precision)) +
      geom_path() +
      theme_bw() +
      facet_grid(~factor(
        Calculator, 
        levels = c("Clinical features", "Clinical features and biomarkers"), 
        labels = c("Early-insulin-treated: clinical features, no biomarkers", 
                   "Early-insulin-treated: clinical features, with biomarkers")), 
        scales = "free",
        labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
      scale_y_continuous("Precision", 
                         labels = scales::percent, 
                         limits = c(0, 1)) +
      scale_x_continuous("Recall", 
                         labels = scales::percent, 
                         limits = c(0, 1)) +
      theme_bw() +
      geom_label(
        data = dat_text %>%
          filter(Dataset == "UNITED" & Model == "Type 1"),
        mapping = aes(x = 0.55, y = 0.91, label = auc_full, hjust = "center"),
        size = 7,
        label.r = unit(0, "pt"),
        label.padding=unit(0.4, "lines")
      ) +
      theme(
        panel.spacing.x = unit(1.5, "lines")
      )
  ),
  
  patchwork::free(
    
    prec_recal_curves %>%
      filter(Dataset == "UNITED" & Model == "Type 2") %>%
      mutate(Calculator = factor(
        Calculator, 
        levels = "No Biomarkers", 
        labels = "Not-early-insulin-treated")) %>%
      ggplot(aes(x = recall, y = precision)) +
      geom_path() +
      theme_bw() +
      facet_grid(~~factor(
        Calculator, 
        levels = c("Not-early-insulin-treated"), 
        labels = c("Not-early-insulin-treated: clinical features, no biomarkers")), 
        scales = "free",
        labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
      scale_y_continuous("Precision", 
                         labels = scales::percent, 
                         limits = c(0, 1)) +
      scale_x_continuous("Recall", 
                         labels = scales::percent, 
                         limits = c(0, 1)) +
      theme_bw() +
      geom_label(
        data = dat_text %>%
          filter(Dataset == "UNITED" & Model == "Type 2") %>%
          mutate(
            Calculator = factor(
              Calculator, levels = c("Clinical features"), 
              labels = c("Not-early-insulin-treated"))
          ),
        mapping = aes(x = 0.55, y = 0.91, label = auc_full, hjust = "center"),
        size = 7,
        label.r = unit(0, "pt"),
        label.padding=unit(0.4, "lines")
      )
    
  ),
  
  nrow = 2, ncol = 1
  
  
) + patchwork::plot_annotation(tag_levels = list(c("A", "B"))) +
  patchwork::plot_layout(
    design = "
    AAAAA
    #BBB#
    "
  ) &
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 11)
  )



pdf("figures/SF_12_UNITED_precision_recall.pdf", width = 10, height = 8)
plot_prob_prec_recal_united
dev.off()

#:------------------------------------------------------------


