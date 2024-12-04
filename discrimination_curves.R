#:--------------------------------------------------------
#   
# In this file we check AUROC for the MODY models
#
#:--------------------------------------------------------

# load libraries
library(tidyverse)
library(nimble)
library(pROC)
library(PRROC)


# load functions needed
source("data/create_data.R")
source("new_data_predictions/prediction_functions.R")

# load files required
predictions_dataset.UNITED_type1_all_genes_no_T_full <- readRDS("model_predictions/predictions_dataset.UNITED_type1_all_genes_no_T_full.rds")
predictions_dataset.UNITED_type1_all_genes_with_T_full <- readRDS("model_predictions/predictions_dataset.UNITED_type1_all_genes_with_T_full.rds")
predictions_dataset.UNITED_type2_all_genes_new_full <- readRDS("model_predictions/predictions_dataset.UNITED_type2_all_genes_new_full.rds")

predictions_dataset.UNITED_type1_all_genes_no_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_all_genes_no_T.rds")
predictions_dataset.UNITED_type1_all_genes_with_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_all_genes_with_T.rds")
predictions_dataset.UNITED_type2_all_genes_new <- readRDS("model_predictions/predictions_dataset.UNITED_type2_all_genes_new.rds")

# load datasets
## Load population representative dataset
dataset.UNITED_type1_all_genes <- create_data(dataset = "united t1d", commonmody = FALSE) %>%
  
  ## if MODY testing missing, change to 0
  mutate(M = ifelse(is.na(M), 0, M))

dataset.UNITED_type2_all_genes <- create_data(dataset = "united t2d", commonmody = FALSE)


#:------------------------------------------------------------

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

### No biomarker models
auc_T1D_no_T_united_all_genes <- calc_auroc(dataset.UNITED_type1_all_genes$M, predictions_dataset.UNITED_type1_all_genes_no_T_full, thinning = 10)
# quantile(auc_T1D_no_T_united_all_genes, probs = c(0.025, 0.5, 0.975)) # thinning = 10
# 2.5%       50%     97.5% 
# 0.7053243 0.7924399 0.8190743 

### Biomarker models
auc_T1D_with_T_united_all_genes <- calc_auroc(dataset.UNITED_type1_all_genes$M, predictions_dataset.UNITED_type1_all_genes_with_T_full, thinning = 10)
# quantile(auc_T1D_with_T_united_all_genes, probs = c(0.025, 0.5, 0.975)) # thinning = 10
# 2.5%       50%     97.5% 
# 0.9487113 0.9768041 0.9779210 

## Type 2 UNITED
auc_T2D_new_united_all_genes <- calc_auroc(dataset.UNITED_type2_all_genes$M, predictions_dataset.UNITED_type2_all_genes_new_full, thinning = 10)
# quantile(auc_T2D_new_united_all_genes, probs = c(0.025, 0.5, 0.975)) # thinning = 10
# 2.5%       50%     97.5% 
# 0.8465473 0.8618926 0.8768116  

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
pr_auc_T1D_no_T_united_all_genes <- calc_auc_pr(dataset.UNITED_type1_all_genes$M, predictions_dataset.UNITED_type1_all_genes_no_T_full, thinning = 10)
# quantile(pr_auc_T1D_no_T_united_all_genes, probs = c(0.025, 0.5, 0.975)) # thinning = 10
# 2.5%        50%      97.5% 
# 0.01548813 0.02644908 0.03477454 

### Biomarker models
pr_auc_T1D_with_T_united_all_genes <- calc_auc_pr(dataset.UNITED_type1_all_genes$M, predictions_dataset.UNITED_type1_all_genes_with_T_full, thinning = 10)
# quantile(pr_auc_T1D_with_T_united_all_genes, probs = c(0.025, 0.5, 0.975)) # thinning = 10
# 2.5%        50%      97.5% 
# 0.09617254 0.19216510 0.22411871 

## Type 2 UNITED
pr_auc_T2D_new_united_all_genes <- calc_auc_pr(dataset.UNITED_type2_all_genes$M, predictions_dataset.UNITED_type2_all_genes_new_full, thinning = 10)
# quantile(pr_auc_T2D_new_united_all_genes, probs = c(0.025, 0.5, 0.975)) # thinning = 10
# 2.5%       50%     97.5% 
# 0.5980896 0.6291710 0.6611516 


# #:------------------------------------------------------------
# 
# # Calculate ROC with intervals
# calc_roc <- function(data, predictions, thinning = 100) {
#   
#   output <- NULL
#   
#   # sequence to iterate through
#   sequence_list <- seq(1, nrow(predictions), thinning)
#   
#   for (i in 1:length(sequence_list)) {
#     
#     # print the current iteration
#     if (i %% 1000 == 0) {
#       print(paste(i, "out of", length(sequence_list)))
#     }
#     
#     ## calculate ROC
#     interim <- pROC::roc(response = data, predictor = predictions[sequence_list[i],], levels = c(0,1), direction = "<") %>%
#       magrittr::extract(c(2:3)) %>%
#       as.data.frame() %>%
#       mutate(iteration = paste0(sequence_list[i]))
#     
#     output <- rbind(output, interim)
#     
#   }
#   
#   return(output)
#   
# }
# 
# 
# ## Type 1 UNITED
# 
# ### No biomarker models
# roc_T1D_no_T_united_all_genes <- calc_roc(dataset.UNITED_type1_all_genes$M, predictions_dataset.UNITED_type1_all_genes_no_T_full, thinning = 100)
# 
# # plot for ROC with grey being iterations, black being the ROC for average prediction
# plot_roc_T1D_no_T_united_all_genes <- ggplot() +
#   ## all iterations
#   geom_path(
#     data = roc_T1D_no_T_united_all_genes,
#     aes(x = 1-specificities, y= sensitivities, group = iteration), colour = "grey"
#   ) +
#   ## average predictions
#   geom_path(
#     data = pROC::roc(response = dataset.UNITED_type1_all_genes$M, predictor = colMeans(predictions_dataset.UNITED_type1_all_genes_no_T_full), levels = c(0,1), direction = "<") %>%
#       magrittr::extract(c(2:3)) %>%
#       as.data.frame(),
#     aes(x = 1-specificities, y= sensitivities), colour = "black"
#   )
# 
# ### Biomarker models
# roc_T1D_with_T_united_all_genes <- calc_roc(dataset.UNITED_type1_all_genes$M, predictions_dataset.UNITED_type1_all_genes_with_T_full, thinning = 100)
# 
# # plot for ROC with grey being iterations, black being the ROC for average prediction
# plot_roc_T1D_with_T_united_all_genes <- ggplot() +
#   ## all iterations
#   geom_path(
#     data = roc_T1D_with_T_united_all_genes,
#     aes(x = 1-specificities, y= sensitivities, group = iteration), colour = "grey"
#   ) +
#   ## average predictions
#   geom_path(
#     data = pROC::roc(response = dataset.UNITED_type1_all_genes$M, predictor = colMeans(predictions_dataset.UNITED_type1_all_genes_with_T_full)) %>%
#       magrittr::extract(c(2:3)) %>%
#       as.data.frame(),
#     aes(x = 1-specificities, y= sensitivities), colour = "black"
#   )
# 
# ## Type 2 UNITED
# roc_T2D_new_united_all_genes <- calc_roc(dataset.UNITED_type2_all_genes$M, predictions_dataset.UNITED_type2_all_genes_new_full, thinning = 100)
# 
# # plot for ROC with grey being iterations, black being the ROC for average prediction
# plot_roc_T2D_new_united_all_genes <- ggplot() +
#   ## all iterations
#   geom_path(
#     data = roc_T2D_new_united_all_genes,
#     aes(x = 1-specificities, y= sensitivities, group = iteration), colour = "grey"
#   ) +
#   ## average predictions
#   geom_path(
#     data = pROC::roc(response = dataset.UNITED_type2_all_genes$M, predictor = colMeans(predictions_dataset.UNITED_type2_all_genes_new_full)) %>%
#       magrittr::extract(c(2:3)) %>%
#       as.data.frame(),
#     aes(x = 1-specificities, y= sensitivities), colour = "black"
#   )
# 
# 
# 
# #:------------------------------------------------------------
# 
# 
# # Calculate Precision-recall with intervals
# calc_prec_recal_curve <- function(data, predictions, thinning = 100) {
#   
#   output <- NULL
#   
#   # sequence to iterate through
#   sequence_list <- seq(1, nrow(predictions), thinning)
#   
#   for (i in 1:length(sequence_list)) {
#     
#     # print the current iteration
#     if (i %% 100 == 0) {
#       print(paste(i, "out of", length(sequence_list)))
#     }
#     
#     ## calculate ROC
#     interim <- pROC::coords(pROC::roc(response = data, predictor = predictions[sequence_list[i],], levels = c(0,1), direction = "<"), ret = c("precision", "recall")) %>%
#       as.data.frame() %>%
#       mutate(iteration = paste0(sequence_list[i]))
#     
#     output <- rbind(output, interim)
#     
#   }
#   
#   return(output)
#   
# }
# 
# 
# ## Type 1 UNITED
# 
# ### No biomarker models
# prec_recal_T1D_no_T_united_all_genes <- calc_prec_recal_curve(dataset.UNITED_type1_all_genes$M, predictions_dataset.UNITED_type1_all_genes_no_T_full, thinning = 100)
# 
# # plot for ROC with grey being iterations, black being the ROC for average prediction
# plot_prec_recal_T1D_no_T_united_all_genes <- ggplot() +
#   ## all iterations
#   geom_path(
#     data = prec_recal_T1D_no_T_united_all_genes,
#     aes(x = recall, y= precision, group = iteration), colour = "grey"
#   ) +
#   ## average predictions
#   geom_path(
#     data = pROC::coords(pROC::roc(response = dataset.UNITED_type1_all_genes$M, predictor = colMeans(predictions_dataset.UNITED_type1_all_genes_no_T_full)), ret = c("precision", "recall")) %>%
#       as.data.frame(),
#     aes(x = recall, y= precision), colour = "black"
#   )
# 
# ### Biomarker models
# prec_recal_T1D_with_T_united_all_genes <- calc_prec_recal_curve(dataset.UNITED_type1_all_genes$M, predictions_dataset.UNITED_type1_all_genes_with_T_full, thinning = 100)
# 
# # plot for ROC with grey being iterations, black being the ROC for average prediction
# plot_prec_recal_T1D_with_T_united_all_genes <- ggplot() +
#   ## all iterations
#   geom_path(
#     data = prec_recal_T1D_with_T_united_all_genes,
#     aes(x = recall, y= precision, group = iteration), colour = "grey"
#   ) +
#   ## average predictions
#   geom_path(
#     data = pROC::coords(pROC::roc(response = dataset.UNITED_type1_all_genes$M, predictor = colMeans(predictions_dataset.UNITED_type1_all_genes_with_T_full)), ret = c("precision", "recall")) %>%
#       as.data.frame(),
#     aes(x = recall, y= precision), colour = "black"
#   )
# 
# ## Type 2 UNITED
# prec_recal_T2D_new_united_all_genes <- calc_prec_recal_curve(dataset.UNITED_type2_all_genes$M, predictions_dataset.UNITED_type2_all_genes_new_full, thinning = 100)
# 
# # plot for ROC with grey being iterations, black being the ROC for average prediction
# plot_prec_recal_T2D_new_united_all_genes <- ggplot() +
#   ## all iterations
#   geom_path(
#     data = prec_recal_T2D_new_united_all_genes,
#     aes(x = recall, y= precision, group = iteration), colour = "grey"
#   ) +
#   ## average predictions
#   geom_path(
#     data = pROC::coords(pROC::roc(response = dataset.UNITED_type2_all_genes$M, predictor = colMeans(predictions_dataset.UNITED_type2_all_genes_new_full)), ret = c("precision", "recall")) %>%
#       as.data.frame(),
#     aes(x = recall, y= precision), colour = "black"
#   )



#:------------------------------------------------------------------------------

# Boxplot and roc curves

roc_curves <- data.frame(prob = colMeans(predictions_dataset.UNITED_type1_all_genes_with_T_full)) %>%
  cbind(Mody = dataset.UNITED_type1_all_genes$M) %>%
  pROC::roc(response = Mody, predictor = prob) %>%
  magrittr::extract(2:3) %>%
  as.data.frame() %>%
  mutate(
    auc =  unname(data.frame(prob = colMeans(predictions_dataset.UNITED_type1_all_genes_with_T_full)) %>%
                       cbind(Mody = dataset.UNITED_type1_all_genes$M) %>%
                       pROC::roc(response = Mody, predictor = prob) %>%
                       magrittr::extract(c(9)) %>%
                       unlist()),
    auc_low = quantile(auc_T1D_with_T_united_all_genes, probs = c(0.025)),
    auc_high = quantile(auc_T1D_with_T_united_all_genes, probs= c(0.975)),
    mean = mean(colMeans(predictions_dataset.UNITED_type1_all_genes_with_T_full), na.rm = TRUE)
  ) %>%
  mutate(Dataset = "UNITED", Model = "Type 1", Calculator = "Biomarkers") %>%
  rbind(
    data.frame(prob = colMeans(predictions_dataset.UNITED_type1_all_genes_no_T_full)) %>%
      cbind(Mody = dataset.UNITED_type1_all_genes$M) %>%
      pROC::roc(response = Mody, predictor = prob) %>%
      magrittr::extract(2:3) %>%
      as.data.frame() %>%
      mutate(
        auc = unname(data.frame(prob = colMeans(predictions_dataset.UNITED_type1_all_genes_no_T_full)) %>%
                          cbind(Mody = dataset.UNITED_type1_all_genes$M) %>%
                          pROC::roc(response = Mody, predictor = prob) %>%
                          magrittr::extract(c(9)) %>%
                          unlist()),
        auc_low = quantile(auc_T1D_no_T_united_all_genes, probs = c(0.025)),
        auc_high = quantile(auc_T1D_no_T_united_all_genes, probs= c(0.975)),
        mean = mean(colMeans(predictions_dataset.UNITED_type1_all_genes_no_T_full), na.rm = TRUE)
      ) %>%
      mutate(Dataset = "UNITED", Model = "Type 1", Calculator = "No Biomarkers"), 
    data.frame(prob = colMeans(predictions_dataset.UNITED_type2_all_genes_new_full)) %>%
      cbind(Mody = dataset.UNITED_type2_all_genes$M) %>%
      pROC::roc(response = Mody, predictor = prob) %>%
      magrittr::extract(2:3) %>%
      as.data.frame() %>%
      mutate(
        auc = unname(data.frame(prob = colMeans(predictions_dataset.UNITED_type2_all_genes_new_full)) %>%
                          cbind(Mody = dataset.UNITED_type2_all_genes$M) %>%
                          pROC::roc(response = Mody, predictor = prob) %>%
                          magrittr::extract(c(9)) %>%
                          unlist()),
        auc_low = quantile(auc_T2D_new_united_all_genes, probs = c(0.025)),
        auc_high = quantile(auc_T2D_new_united_all_genes, probs= c(0.975)),
        mean = mean(colMeans(predictions_dataset.UNITED_type2_all_genes_new_full), na.rm = TRUE)
      ) %>%
      mutate(Dataset = "UNITED", Model = "Type 2", Calculator = "No Biomarkers")
  )


dat_text <- roc_curves %>%
  select(-sensitivities, -specificities) %>%
  distinct() %>%
  mutate(
    auc_full = paste0("AUC: ", signif(auc, 2), " [", signif(auc_low, 2), "-", signif(auc_high, 2), "]"),
    mean = paste0("Mean prob:", signif(mean, 2)*100, "%"),
    Calculator = factor(Calculator, levels = c("Biomarkers", "No Biomarkers"), labels = c("Clinical features and biomarkers", "Clinical features"))
  )


plot_prob_density_rocs_external_bio_only <- patchwork::wrap_plots(
  # Panel A - External insulin-treated
  patchwork::wrap_plots(
    patchwork::wrap_plots(
      # Density
      dataset.UNITED_type1_all_genes %>%
        mutate(prob = predictions_dataset.UNITED_type1_all_genes_with_T$prob) %>%
        filter(M == 0) %>%
        select(M, prob) %>%
        rename("Mody" = "M") %>%
        mutate(
          Mody = factor(Mody, levels = c(0, 1), labels = c("Non-MODY", "MODY")),
        ) %>%
        ggplot() +
        geom_density(aes(x = prob), fill = "grey") +
        coord_cartesian(xlim =c(0, 1)) +
        scale_x_continuous(labels = scales::percent) +
        theme_bw() +
        theme(
          panel.border = element_blank(),
          axis.line = element_line(),
          axis.ticks.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.title.x = element_blank()
        ) +
        ylab("Non-MODY (n=935)"),
      #point
      dataset.UNITED_type1_all_genes %>%
        mutate(prob = predictions_dataset.UNITED_type1_all_genes_with_T$prob) %>%
        filter(M == 1) %>%
        select(M, prob) %>%
        rename("Mody" = "M") %>%
        mutate(
          Mody = factor(Mody, levels = c(0, 1), labels = c("Non-MODY", "MODY")),
        ) %>%
        ggplot() +
        geom_point(aes(x = prob, y=0), position = position_jitter(height = 0.1, seed = 20)) +
        coord_cartesian(xlim =c(0, 1), ylim = c(-0.15, 0.15)) +
        scale_x_continuous(labels = scales::percent) +
        theme_classic() +
        theme(
          axis.ticks.y = element_blank(), 
          axis.text.y = element_blank(), 
          #axis.title.y = element_blank()
        ) +
        ylab("MODY \n cases \n (n=5)") +
        xlab("Model probabilities"),
      ncol = 1, nrow = 2, heights = c(4,1)
    ),
    #rocs
    roc_curves %>%
      filter(Dataset == "UNITED" & Model == "Type 1" & Calculator == "Biomarkers") %>%
      mutate(
        Calculator = factor(Calculator, levels = c("No Biomarkers", "Biomarkers"), labels = c("Clinical features", "Clinical features and biomarkers")),
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
          filter(Dataset == "UNITED" & Model == "Type 1" & Calculator == "Clinical features and biomarkers"),
        mapping = aes(x = 0.55, y = 0.1, label = auc_full, hjust = "center"),
        size = 7,
        label.r = unit(0, "pt"),
        label.padding=unit(0.4, "lines")
      ) +
      theme(
        panel.spacing.x = unit(1.5, "lines")
      ), 
    ncol = 2, nrow = 1
  ),
  # Panel B - External non-insulin-treated
  patchwork::wrap_plots(
    patchwork::wrap_plots(
      # Density
      dataset.UNITED_type2_all_genes %>%
        mutate(prob = predictions_dataset.UNITED_type2_all_genes_new$prob) %>%
        filter(M == 0) %>%
        select(M, prob) %>%
        rename("Mody" = "M") %>%
        mutate(
          Mody = factor(Mody, levels = c(0, 1), labels = c("Non-MODY", "MODY")),
        ) %>%
        ggplot() +
        geom_density(aes(x = prob), fill = "grey") +
        coord_cartesian(xlim =c(0, 1)) +
        scale_x_continuous(labels = scales::percent) +
        theme_bw() +
        theme(
          panel.border = element_blank(),
          axis.line = element_line(),
          axis.ticks.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.title.x = element_blank()
        ) +
        ylab("Non-MODY (n=70)"),
      #point
      dataset.UNITED_type2_all_genes %>%
        mutate(prob = predictions_dataset.UNITED_type2_all_genes_new$prob) %>%
        filter(M == 1) %>%
        select(M, prob) %>%
        rename("Mody" = "M") %>%
        mutate(
          Mody = factor(Mody, levels = c(0, 1), labels = c("Non-MODY", "MODY")),
        ) %>%
        ggplot() +
        geom_point(aes(x=prob, y=0), position = position_jitter(height = 0.1, seed = 20)) +
        coord_cartesian(xlim =c(0, 1), ylim = c(-0.15, 0.15)) +
        scale_x_continuous(labels = scales::percent) +
        theme_classic() +
        theme(
          axis.ticks.y = element_blank(), 
          axis.text.y = element_blank(), 
          #axis.title.y = element_blank()
        ) +
        ylab("MODY \n cases \n (n=15)") +
        xlab("Model probabilities"),
      ncol = 1, nrow = 2, heights = c(4,1)
    ),
    #rocs
    
    roc_curves %>%
      filter(Dataset == "UNITED" & Model == "Type 2") %>%
      mutate(Calculator = factor(Calculator, levels = "No Biomarkers", labels = "Not-early-insulin-treated")) %>%
      ggplot(aes(x = 1- specificities, y = sensitivities)) +
      geom_path() +
      theme_bw() +
      scale_y_continuous("Sensitivity", labels = scales::percent) +
      scale_x_continuous("1- Specificity", labels = scales::percent) +
      theme_bw() +
      geom_label(
        data = dat_text %>%
          filter(Dataset == "UNITED" & Model == "Type 2") %>%
          mutate(
            Calculator = factor(Calculator, levels = c("Clinical features"), labels = c("Not-early-insulin-treated"))
          ),
        mapping = aes(x = 0.55, y = 0.1, label = auc_full, hjust = "center"),
        size = 7,
        label.r = unit(0, "pt"),
        label.padding=unit(0.4, "lines")
      ),
    
    ncol = 2, nrow = 1
  ),
  
  ncol = 1
  
) + patchwork::plot_annotation(tag_levels = list(c("A.1", "", "A.2", "B.1", "","B.2"))) &
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 11)
  )

# Making plots
pdf("figures/united_boxplot_roc.pdf", width = 13, height = 9)
plot_prob_density_rocs_external_bio_only
dev.off()


#:-------------------------------------------------------------


plot_prob_rocs_united <- patchwork::wrap_plots(
  
  # ROC T1D
  patchwork::free(
    
    roc_curves %>%
      filter(Dataset == "UNITED" & Model == "Type 1") %>%
      mutate(
        Calculator = factor(Calculator, levels = c("No Biomarkers", "Biomarkers"), labels = c("Clinical features", "Clinical features and biomarkers")),
        iteration = 0
      ) %>%
      ggplot(aes(x = 1- specificities, y = sensitivities)) +
      geom_path() +
      theme_bw() +
      facet_grid(~factor(Calculator, levels = c("Clinical features", "Clinical features and biomarkers"), labels = c("Early-insulin-treated: clinical features, no biomarkers", "Early-insulin-treated: clinical features, with biomarkers")), scales = "free",labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
      scale_y_continuous("Sensitivity", labels = scales::percent) +
      scale_x_continuous("1- Specificity", labels = scales::percent) +
      theme_bw() +
      geom_label(
        data = dat_text %>%
          filter(Dataset == "UNITED" & Model == "Type 1"),
        mapping = aes(x = 0.55, y = 0.1, label = auc_full, hjust = "center"),
        size = 7,
        label.r = unit(0, "pt"),
        label.padding=unit(0.4, "lines")
      ) +
      theme(
        panel.spacing.x = unit(1.5, "lines")
      )
  ),
  
  patchwork::free(
    
    roc_curves %>%
      filter(Dataset == "UNITED" & Model == "Type 2") %>%
      mutate(Calculator = factor(Calculator, levels = "No Biomarkers", labels = "Not-early-insulin-treated")) %>%
      ggplot(aes(x = 1- specificities, y = sensitivities)) +
      geom_path() +
      theme_bw() +
      facet_grid(~~factor(Calculator, levels = c("Not-early-insulin-treated"), labels = c("Not-early-insulin-treated: clinical features, no biomarkers")), scales = "free",labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
      scale_y_continuous("Sensitivity", labels = scales::percent) +
      scale_x_continuous("1- Specificity", labels = scales::percent) +
      theme_bw() +
      geom_label(
        data = dat_text %>%
          filter(Dataset == "UNITED" & Model == "Type 2") %>%
          mutate(
            Calculator = factor(Calculator, levels = c("Clinical features"), labels = c("Not-early-insulin-treated"))
          ),
        mapping = aes(x = 0.55, y = 0.1, label = auc_full, hjust = "center"),
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


pdf("figures/united_roc.pdf", width = 10, height = 8)
plot_prob_rocs_united
dev.off()



################


prec_recal_curves <- data.frame(prob = colMeans(predictions_dataset.UNITED_type1_all_genes_with_T_full)) %>%
  cbind(Mody = dataset.UNITED_type1_all_genes$M) %>%
  pROC::roc(response = Mody, predictor = prob) %>%
  pROC::coords(ret = c("precision", "recall")) %>%
  as.data.frame() %>%
  mutate(
    auc = calc_auc_pr(dataset.UNITED_type1_all_genes$M, t(as.data.frame(colMeans(predictions_dataset.UNITED_type1_all_genes_with_T_full))), thinning = 1),
    auc_low = quantile(pr_auc_T1D_with_T_united_all_genes, probs = c(0.025)),
    auc_high = quantile(pr_auc_T1D_with_T_united_all_genes, probs= c(0.975)),
    mean = mean(colMeans(predictions_dataset.UNITED_type1_all_genes_with_T_full), na.rm = TRUE)
  ) %>%
  mutate(Dataset = "UNITED", Model = "Type 1", Calculator = "Biomarkers") %>%
  rbind(
    data.frame(prob = colMeans(predictions_dataset.UNITED_type1_all_genes_no_T_full)) %>%
      cbind(Mody = dataset.UNITED_type1_all_genes$M) %>%
      pROC::roc(response = Mody, predictor = prob) %>%
      pROC::coords(ret = c("precision", "recall")) %>%
      as.data.frame() %>%
      mutate(
        auc = calc_auc_pr(dataset.UNITED_type1_all_genes$M, t(as.data.frame(colMeans(predictions_dataset.UNITED_type1_all_genes_no_T_full))), thinning = 1),
        auc_low = quantile(pr_auc_T1D_no_T_united_all_genes, probs = c(0.025)),
        auc_high = quantile(pr_auc_T1D_no_T_united_all_genes, probs= c(0.975)),
        mean = mean(colMeans(predictions_dataset.UNITED_type1_all_genes_no_T_full), na.rm = TRUE)
      ) %>%
      mutate(Dataset = "UNITED", Model = "Type 1", Calculator = "No Biomarkers"), 
    data.frame(prob = colMeans(predictions_dataset.UNITED_type2_all_genes_new_full)) %>%
      cbind(Mody = dataset.UNITED_type2_all_genes$M) %>%
      pROC::roc(response = Mody, predictor = prob) %>%
      pROC::coords(ret = c("precision", "recall")) %>%
      as.data.frame() %>%
      mutate(
        auc = calc_auc_pr(dataset.UNITED_type2_all_genes$M, t(as.data.frame(colMeans(predictions_dataset.UNITED_type2_all_genes_new_full))), thinning = 1),
        auc_low = quantile(pr_auc_T2D_new_united_all_genes, probs = c(0.025)),
        auc_high = quantile(pr_auc_T2D_new_united_all_genes, probs= c(0.975)),
        mean = mean(colMeans(predictions_dataset.UNITED_type2_all_genes_new_full), na.rm = TRUE)
      ) %>%
      mutate(Dataset = "UNITED", Model = "Type 2", Calculator = "No Biomarkers")
  )


dat_text <- prec_recal_curves %>%
  select(-precision, -recall) %>%
  distinct() %>%
  mutate(
    auc_full = paste0("AUC: ", signif(auc, 2), " [", signif(auc_low, 2), "-", signif(auc_high, 2), "]"),
    mean = paste0("Mean prob:", signif(mean, 2)*100, "%"),
    Calculator = factor(Calculator, levels = c("Biomarkers", "No Biomarkers"), labels = c("Clinical features and biomarkers", "Clinical features"))
  )



plot_prob_prec_recal_united <- patchwork::wrap_plots(
  
  # ROC T1D
  patchwork::free(
    
    prec_recal_curves %>%
      filter(Dataset == "UNITED" & Model == "Type 1") %>%
      mutate(
        Calculator = factor(Calculator, levels = c("No Biomarkers", "Biomarkers"), labels = c("Clinical features", "Clinical features and biomarkers")),
        iteration = 0
      ) %>%
      ggplot(aes(x = recall, y = precision)) +
      geom_path() +
      theme_bw() +
      facet_grid(~factor(Calculator, levels = c("Clinical features", "Clinical features and biomarkers"), labels = c("Early-insulin-treated: clinical features, no biomarkers", "Early-insulin-treated: clinical features, with biomarkers")), scales = "free",labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
      scale_y_continuous("Precision", labels = scales::percent, limits = c(0, 1)) +
      scale_x_continuous("Recall", labels = scales::percent, limits = c(0, 1)) +
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
      mutate(Calculator = factor(Calculator, levels = "No Biomarkers", labels = "Not-early-insulin-treated")) %>%
      ggplot(aes(x = recall, y = precision)) +
      geom_path() +
      theme_bw() +
      facet_grid(~~factor(Calculator, levels = c("Not-early-insulin-treated"), labels = c("Not-early-insulin-treated: clinical features, no biomarkers")), scales = "free",labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
      scale_y_continuous("Precision", labels = scales::percent, limits = c(0, 1)) +
      scale_x_continuous("Recall", labels = scales::percent, limits = c(0, 1)) +
      theme_bw() +
      geom_label(
        data = dat_text %>%
          filter(Dataset == "UNITED" & Model == "Type 2") %>%
          mutate(
            Calculator = factor(Calculator, levels = c("Clinical features"), labels = c("Not-early-insulin-treated"))
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



pdf("figures/united_prec_recal_thin.pdf", width = 10, height = 8)
plot_prob_prec_recal_united
dev.off()









































































# plot_prob_rocs_united <- patchwork::wrap_plots(
#   
#   # roc insulin-treated
#   patchwork::free(roc_curves %>%
#                     filter(Dataset == "UNITED" & Model == "Type 1") %>%
#                     mutate(
#                       Calculator = factor(Calculator, levels = c("No Biomarkers", "Biomarkers"), labels = c("Clinical features", "Clinical features and biomarkers")),
#                       iteration = 0
#                     ) %>%
#                     ggplot(aes(x = 1- specificities, y = sensitivities)) +
#                     geom_path(
#                       data = roc_T1D_no_T_united_all_genes %>%
#                         mutate(
#                           Calculator = "Clinical features"
#                         ) %>%
#                         rbind(
#                           roc_T1D_with_T_united_all_genes %>%
#                             mutate(Calculator = "Clinical features and biomarkers")
#                         ),
#                       aes(group = iteration), colour = "grey"
#                     ) + 
#                     geom_path() +
#                     theme_bw() +
#                     facet_grid(~factor(Calculator, levels = c("Clinical features", "Clinical features and biomarkers"), labels = c("Early-insulin-treated: clinical features, no biomarkers", "Early-insulin-treated: clinical features, with biomarkers")), scales = "free",labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
#                     scale_y_continuous("Sensitivity", labels = scales::percent) +
#                     scale_x_continuous("1- Specificity", labels = scales::percent) +
#                     theme_bw() +
#                     geom_label(
#                       data = roc_curves %>%
#                         select(-sensitivities, -specificities) %>%
#                         distinct() %>%
#                         mutate(
#                           auc = paste0(" AUC:", signif(auc, 2), " "),
#                           mean = paste0("Mean prob:", signif(mean, 2)*100, "%")
#                         ) %>%
#                         filter(Dataset == "UNITED" & Model == "Type 1") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = c("No Biomarkers", "Biomarkers"), labels = c("Clinical features", "Clinical features and biomarkers"))
#                         ),
#                       mapping = aes(x = -Inf, y = -Inf, label = auc),
#                       size = 7,
#                       label.size = NA,
#                       hjust = -0.8,
#                       vjust = -0.5
#                     ) +
#                     theme(
#                       panel.spacing.x = unit(1.5, "lines")
#                     )),
#   # roc not early insulin treated
#   patchwork::free(roc_curves %>%
#                     filter(Dataset == "UNITED" & Model == "Type 2") %>%
#                     mutate(Calculator = factor(Calculator, levels = "No Biomarkers", labels = "Not-early-insulin-treated")) %>%
#                     ggplot(aes(x = 1- specificities, y = sensitivities)) +
#                     geom_path(
#                       data = roc_T2D_new_united_all_genes %>%
#                         mutate(
#                           Calculator = factor("Not-early-insulin-treated")
#                         ),
#                       aes(group = iteration), colour = "grey"
#                     ) +
#                     geom_path() +
#                     theme_bw() +
#                     facet_grid(~factor(Calculator, levels = "Not-early-insulin-treated", labels = "Not-early-insulin-treated: clinical features, no biomarkers"), scales = "free",labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
#                     scale_y_continuous("Sensitivity", labels = scales::percent) +
#                     scale_x_continuous("1- Specificity", labels = scales::percent) +
#                     theme_bw() +
#                     geom_label(
#                       data = dat_text %>%
#                         filter(Dataset == "UNITED" & Model == "Type 2") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = c("Clinical features"), labels = c("Not-early-insulin-treated"))
#                         ),
#                       mapping = aes(x = -Inf, y = -Inf, label = auc),
#                       size = 7,
#                       label.size = NA,
#                       hjust = -0.8,
#                       vjust = -0.5
#                     )),
#   
#   ncol = 1
#   
# ) + patchwork::plot_annotation(tag_levels = list(c("A", "B"))) &
#   theme(
#     axis.text = element_text(size = 14),
#     axis.title = element_text(size = 16),
#     strip.text = element_text(size = 11)
#   )
# 
# 
# plot_prob_rocs_united1 <- patchwork::wrap_plots(
#   
#   # roc insulin-treated
#   patchwork::free(roc_curves %>%
#                     filter(Dataset == "UNITED" & Model == "Type 1") %>%
#                     mutate(
#                       Calculator = factor(Calculator, levels = c("No Biomarkers", "Biomarkers"), labels = c("Clinical features", "Clinical features and biomarkers")),
#                       iteration = 0
#                     ) %>%
#                     ggplot(aes(x = 1- specificities, y = sensitivities)) +
#                     geom_path(
#                       data = roc_T1D_no_T_united_all_genes %>%
#                         mutate(
#                           Calculator = "Clinical features"
#                         ) %>%
#                         rbind(
#                           roc_T1D_with_T_united_all_genes %>%
#                             mutate(Calculator = "Clinical features and biomarkers")
#                         ),
#                       aes(group = iteration), colour = "grey"
#                     ) + 
#                     geom_path() +
#                     theme_bw() +
#                     facet_grid(~factor(Calculator, levels = c("Clinical features", "Clinical features and biomarkers"), labels = c("Early-insulin-treated: Clinical features", "Early-insulin-treated: Clinical features and biomarkers")), scales = "free",labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
#                     scale_y_continuous("Sensitivity", labels = scales::percent) +
#                     scale_x_continuous("1- Specificity", labels = scales::percent) +
#                     theme_bw() +
#                     geom_label(
#                       data = roc_curves %>%
#                         select(-sensitivities, -specificities) %>%
#                         distinct() %>%
#                         mutate(
#                           auc = paste0(" AUC:", signif(auc, 2), " "),
#                           mean = paste0("Mean prob:", signif(mean, 2)*100, "%")
#                         ) %>%
#                         filter(Dataset == "UNITED" & Model == "Type 1") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = c("Biomarkers", "No Biomarkers"), labels = c("Clinical features and biomarkers", "Clinical features"))
#                         ),
#                       mapping = aes(x = -Inf, y = -Inf, label = auc),
#                       size = 7,
#                       label.size = NA,
#                       hjust = -0.8,
#                       vjust = -0.5
#                     ) +
#                     theme(
#                       panel.spacing.x = unit(1.5, "lines")
#                     )),
#   # roc not early insulin treated
#   patchwork::free(roc_curves %>%
#                     filter(Dataset == "UNITED" & Model == "Type 2") %>%
#                     mutate(Calculator = factor(Calculator, levels = "No Biomarkers", labels = "Not-early-insulin-treated: Clinical features")) %>%
#                     ggplot(aes(x = 1- specificities, y = sensitivities)) +
#                     geom_path(
#                       data = roc_T2D_new_united_all_genes %>%
#                         mutate(
#                           Calculator = factor("Not-early-insulin-treated: Clinical features")
#                         ),
#                       aes(group = iteration), colour = "grey"
#                     ) +
#                     geom_path() +
#                     theme_bw() +
#                     facet_grid(~Calculator, labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
#                     scale_y_continuous("Sensitivity", labels = scales::percent) +
#                     scale_x_continuous("1- Specificity", labels = scales::percent) +
#                     theme_bw() +
#                     geom_label(
#                       data = dat_text %>%
#                         filter(Dataset == "UNITED" & Model == "Type 2") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = c("Clinical features"), labels = c("Not-early-insulin-treated: Clinical features"))
#                         ),
#                       mapping = aes(x = -Inf, y = -Inf, label = auc),
#                       size = 7,
#                       label.size = NA,
#                       hjust = -0.8,
#                       vjust = -0.5
#                     )),
#   
#   ncol = 1
#   
# ) + patchwork::plot_annotation(tag_levels = list(c("A", "B"))) &
#   theme(
#     axis.text = element_text(size = 14),
#     axis.title = element_text(size = 16),
#     strip.text = element_text(size = 11)
#   )
# 
# plot_prob_rocs_united2 <- patchwork::wrap_plots(
#   
#   # roc insulin-treated
#   patchwork::free(roc_curves %>%
#                     filter(Dataset == "UNITED" & Model == "Type 1") %>%
#                     mutate(
#                       Calculator = factor(Calculator, levels = c("No Biomarkers", "Biomarkers"), labels = c("Clinical features", "Clinical features and biomarkers")),
#                       iteration = 0
#                     ) %>%
#                     ggplot(aes(x = 1- specificities, y = sensitivities)) +
#                     geom_path(
#                       data = roc_T1D_no_T_united_all_genes %>%
#                         mutate(
#                           Calculator = "Clinical features"
#                         ) %>%
#                         rbind(
#                           roc_T1D_with_T_united_all_genes %>%
#                             mutate(Calculator = "Clinical features and biomarkers")
#                         ),
#                       aes(group = iteration), colour = "grey"
#                     ) + 
#                     geom_path() +
#                     theme_bw() +
#                     facet_grid(~factor(Calculator, levels = c("Clinical features", "Clinical features and biomarkers"), labels = c("Early-insulin-treated: Clinical features", "Early-insulin-treated: Clinical features and biomarkers")), scales = "free",labeller = label_wrap_gen(width = 40, multi_line =  TRUE)) +
#                     scale_y_continuous("Sensitivity", labels = scales::percent) +
#                     scale_x_continuous("1- Specificity", labels = scales::percent) +
#                     theme_bw() +
#                     geom_label(
#                       data = roc_curves %>%
#                         select(-sensitivities, -specificities) %>%
#                         distinct() %>%
#                         mutate(
#                           auc = paste0(" AUC:", signif(auc, 2), " "),
#                           mean = paste0("Mean prob:", signif(mean, 2)*100, "%")
#                         ) %>%
#                         filter(Dataset == "UNITED" & Model == "Type 1") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = c("Biomarkers", "No Biomarkers"), labels = c("Clinical features and biomarkers", "Clinical features"))
#                         ),
#                       mapping = aes(x = -Inf, y = -Inf, label = auc),
#                       size = 7,
#                       label.size = NA,
#                       hjust = -0.8,
#                       vjust = -0.5
#                     ) +
#                     theme(
#                       panel.spacing.x = unit(1.5, "lines")
#                     )),
#   # roc not early insulin treated
#   patchwork::free(roc_curves %>%
#                     filter(Dataset == "UNITED" & Model == "Type 2") %>%
#                     mutate(Calculator = factor(Calculator, levels = "No Biomarkers", labels = "Not-early-insulin-treated")) %>%
#                     ggplot(aes(x = 1- specificities, y = sensitivities)) +
#                     geom_path(
#                       data = roc_T2D_new_united_all_genes %>%
#                         mutate(
#                           Calculator = factor("Not-early-insulin-treated")
#                         ),
#                       aes(group = iteration), colour = "grey"
#                     ) +
#                     geom_path() +
#                     theme_bw() +
#                     facet_grid(~Calculator) +
#                     scale_y_continuous("Sensitivity", labels = scales::percent) +
#                     scale_x_continuous("1- Specificity", labels = scales::percent) +
#                     theme_bw() +
#                     geom_label(
#                       data = dat_text %>%
#                         filter(Dataset == "UNITED" & Model == "Type 2") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = c("Clinical features"), labels = c("Not-early-insulin-treated"))
#                         ),
#                       mapping = aes(x = -Inf, y = -Inf, label = auc),
#                       size = 7,
#                       label.size = NA,
#                       hjust = -0.8,
#                       vjust = -0.5
#                     )),
#   
#   ncol = 1
#   
# ) + patchwork::plot_annotation(tag_levels = list(c("A", "B"))) &
#   theme(
#     axis.text = element_text(size = 14),
#     axis.title = element_text(size = 16),
#     strip.text = element_text(size = 11)
#   )
# 
# plot_prob_rocs_united3 <- patchwork::wrap_plots(
#   
#   # roc insulin-treated
#   patchwork::free(roc_curves %>%
#                     filter(Dataset == "UNITED" & Model == "Type 1") %>%
#                     mutate(
#                       Calculator = factor(Calculator, levels = c("No Biomarkers", "Biomarkers"), labels = c("Clinical features", "Clinical features and biomarkers")),
#                       iteration = 0
#                     ) %>%
#                     ggplot(aes(x = 1- specificities, y = sensitivities)) +
#                     geom_path(
#                       data = roc_T1D_no_T_united_all_genes %>%
#                         mutate(
#                           Calculator = "Clinical features"
#                         ) %>%
#                         rbind(
#                           roc_T1D_with_T_united_all_genes %>%
#                             mutate(Calculator = "Clinical features and biomarkers")
#                         ),
#                       aes(group = iteration), colour = "grey"
#                     ) + 
#                     geom_path() +
#                     theme_bw() +
#                     facet_grid(~factor(Calculator, levels = c("Clinical features", "Clinical features and biomarkers"), labels = c("Early-insulin-treated: Clinical features", "Early-insulin-treated: Clinical features and biomarkers")), scales = "free",labeller = label_wrap_gen(width = 40, multi_line =  TRUE)) +
#                     scale_y_continuous("Sensitivity", labels = scales::percent) +
#                     scale_x_continuous("1- Specificity", labels = scales::percent) +
#                     theme_bw() +
#                     geom_label(
#                       data = roc_curves %>%
#                         select(-sensitivities, -specificities) %>%
#                         distinct() %>%
#                         mutate(
#                           auc = paste0(" AUC:", signif(auc, 2), " "),
#                           mean = paste0("Mean prob:", signif(mean, 2)*100, "%")
#                         ) %>%
#                         filter(Dataset == "UNITED" & Model == "Type 1") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = c("Biomarkers", "No Biomarkers"), labels = c("Clinical features and biomarkers", "Clinical features"))
#                         ),
#                       mapping = aes(x = -Inf, y = -Inf, label = auc),
#                       size = 7,
#                       label.size = NA,
#                       hjust = -0.8,
#                       vjust = -0.5
#                     ) +
#                     theme(
#                       panel.spacing.x = unit(1.5, "lines")
#                     )),
#   # roc not early insulin treated
#   patchwork::free(roc_curves %>%
#                     filter(Dataset == "UNITED" & Model == "Type 2") %>%
#                     mutate(Calculator = factor(Calculator, levels = "No Biomarkers", labels = "Not-early-insulin-treated: Clinical features")) %>%
#                     ggplot(aes(x = 1- specificities, y = sensitivities)) +
#                     geom_path(
#                       data = roc_T2D_new_united_all_genes %>%
#                         mutate(
#                           Calculator = factor("Not-early-insulin-treated: Clinical features")
#                         ),
#                       aes(group = iteration), colour = "grey"
#                     ) +
#                     geom_path() +
#                     theme_bw() +
#                     facet_grid(~Calculator, labeller = label_wrap_gen(width = 40, multi_line =  TRUE)) +
#                     scale_y_continuous("Sensitivity", labels = scales::percent) +
#                     scale_x_continuous("1- Specificity", labels = scales::percent) +
#                     theme_bw() +
#                     geom_label(
#                       data = dat_text %>%
#                         filter(Dataset == "UNITED" & Model == "Type 2") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = c("Clinical features"), labels = c("Not-early-insulin-treated: Clinical features"))
#                         ),
#                       mapping = aes(x = -Inf, y = -Inf, label = auc),
#                       size = 7,
#                       label.size = NA,
#                       hjust = -0.8,
#                       vjust = -0.5
#                     )),
#   
#   ncol = 1
#   
# ) + patchwork::plot_annotation(tag_levels = list(c("A", "B"))) &
#   theme(
#     axis.text = element_text(size = 14),
#     axis.title = element_text(size = 16),
#     strip.text = element_text(size = 11)
#   )
# 
# pdf("figures/united_roc_thin_100.pdf", width = 9, height = 8)
# plot_prob_rocs_united +
#   patchwork::plot_layout(
#     design = "
#     AAAA
#     #BB#
#     "
#   )
# plot_prob_rocs_united1 +
#   patchwork::plot_layout(
#     design = "
#     AAAA
#     #BB#
#     "
#   )
# plot_prob_rocs_united2 +
#   patchwork::plot_layout(
#     design = "
#     AAAA
#     #BB#
#     "
#   )
# plot_prob_rocs_united3 +
#   patchwork::plot_layout(
#     design = "
#     AAAA
#     #BB#
#     "
#   )
# dev.off()
# 
# 
# 
# 
# plot_prob_boxplot_united <- patchwork::wrap_plots(
#   
#   # Boxplots
#   patchwork::free(dataset.UNITED_type1_all_genes %>%
#                     select(M) %>%
#                     rename("Mody" = "M") %>%
#                     cbind(
#                       prob_with = colMeans(predictions_dataset.UNITED_type1_all_genes_with_T_full),
#                       prob_without = colMeans(predictions_dataset.UNITED_type1_all_genes_no_T_full)
#                     ) %>%
#                     gather("key", "Probability", -Mody) %>% 
#                     mutate(
#                       Mody = factor(Mody, levels = c(0, 1), labels = c("Negative", "Positive")),
#                       key = factor(key, levels = c("prob_without", "prob_with"), labels = c("Clinical features", "Clinical features and biomarkers"))
#                     ) %>%
#                     ggplot() +
#                     geom_boxplot(aes(y = Probability, x = Mody), colour = c("black", "white", "black", "white"), alpha = c(1, 0, 1, 0)) +
#                     geom_point(aes(y = Probability, x = Mody, colour = Mody, alpha = Mody)) +
#                     scale_y_continuous(labels = scales::percent) +
#                     scale_colour_manual(values = c("white", "black")) +
#                     scale_alpha_manual(values = c(0, 1)) +
#                     facet_wrap(~key) +
#                     theme_bw() +
#                     theme(
#                       legend.position = "none"
#                     )),
#   
#   # Boxplots
#   patchwork::free(dataset.UNITED_type2_all_genes %>%
#                     select(M) %>%
#                     rename("Mody" = "M") %>%
#                     cbind(
#                       Probability = colMeans(predictions_dataset.UNITED_type2_all_genes_new_full),
#                       key = ""
#                     ) %>%
#                     mutate(
#                       Mody = factor(Mody, levels = c(0, 1), labels = c("Negative", "Positive"))
#                     ) %>%
#                     ggplot() +
#                     geom_boxplot(aes(y = Probability, x = Mody)) +
#                     scale_y_continuous(labels = scales::percent) +
#                     theme_bw())
#   
#   , ncol = 1
#   
# ) + patchwork::plot_annotation(tag_levels = list(c("A", "B"))) &
#   theme(
#     axis.text = element_text(size = 14),
#     axis.title = element_text(size = 16),
#     strip.text = element_text(size = 11)
#   )
# 
# 
# pdf("figures/united_boxplot_thin_100.pdf", width = 8, height = 9)
# plot_prob_boxplot_united +
#   patchwork::plot_layout(
#     design = "
#     AAAA
#     #BB#
#     "
#   )
# dev.off()
# 
# 
# ################
# 
# 
# prec_recal_curves <- data.frame(prob = colMeans(predictions_dataset.UNITED_type1_all_genes_with_T_full)) %>%
#   cbind(Mody = dataset.UNITED_type1_all_genes$M) %>%
#   pROC::roc(response = Mody, predictor = prob) %>%
#   pROC::coords(ret = c("precision", "recall")) %>%
#   as.data.frame() %>%
#   mutate(
#     ROCAUC = calc_auc_pr(dataset.UNITED_type1_all_genes$M, t(as.data.frame(colMeans(predictions_dataset.UNITED_type1_all_genes_with_T_full))), thinning = 1),
#     mean = mean(colMeans(predictions_dataset.UNITED_type1_all_genes_with_T_full), na.rm = TRUE)
#   ) %>%
#   mutate(Dataset = "UNITED", Model = "Type 1", Calculator = "Biomarkers") %>%
#   rbind(
#     data.frame(prob = colMeans(predictions_dataset.UNITED_type1_all_genes_no_T_full)) %>%
#       cbind(Mody = dataset.UNITED_type1_all_genes$M) %>%
#       pROC::roc(response = Mody, predictor = prob) %>%
#       pROC::coords(ret = c("precision", "recall")) %>%
#       as.data.frame() %>%
#       mutate(
#         ROCAUC = calc_auc_pr(dataset.UNITED_type1_all_genes$M, t(as.data.frame(colMeans(predictions_dataset.UNITED_type1_all_genes_no_T_full))), thinning = 1),
#         mean = mean(colMeans(predictions_dataset.UNITED_type1_all_genes_no_T_full), na.rm = TRUE)
#       ) %>%
#       mutate(Dataset = "UNITED", Model = "Type 1", Calculator = "No Biomarkers"), 
#     data.frame(prob = colMeans(predictions_dataset.UNITED_type2_all_genes_new_full)) %>%
#       cbind(Mody = dataset.UNITED_type2_all_genes$M) %>%
#       pROC::roc(response = Mody, predictor = prob) %>%
#       pROC::coords(ret = c("precision", "recall")) %>%
#       as.data.frame() %>%
#       mutate(
#         ROCAUC = calc_auc_pr(dataset.UNITED_type2_all_genes$M, t(as.data.frame(colMeans(predictions_dataset.UNITED_type2_all_genes_new_full))), thinning = 1),
#         mean = mean(colMeans(predictions_dataset.UNITED_type2_all_genes_new_full), na.rm = TRUE)
#       ) %>%
#       mutate(Dataset = "UNITED", Model = "Type 2", Calculator = " ")
#   ) %>%
#   rename("auc" = "ROCAUC")
# 
# 
# dat_text <- prec_recal_curves %>%
#   select(-precision, -recall) %>%
#   distinct() %>%
#   mutate(
#     auc = paste0(" AUC:", signif(auc, 2), " "),
#     mean = paste0("Mean prob:", signif(mean, 2)*100, "%"),
#     Calculator = factor(Calculator, levels = c("Biomarkers", "No Biomarkers"), labels = c("Clinical features and biomarkers", "Clinical features"))
#   )
# 
# 
# 
# 
# plot_prob_prec_recal_united <- patchwork::wrap_plots(
#   
#   # roc insulin-treated
#   patchwork::free(prec_recal_curves %>%
#                     filter(Dataset == "UNITED" & Model == "Type 1") %>%
#                     mutate(
#                       Calculator = factor(Calculator, levels = c("No Biomarkers", "Biomarkers"), labels = c("Clinical features", "Clinical features and biomarkers")),
#                       iteration = 0
#                     ) %>%
#                     ggplot(aes(x = recall, y = precision)) +
#                     geom_path(
#                       data = prec_recal_T1D_no_T_united_all_genes %>%
#                         mutate(
#                           Calculator = "Clinical features"
#                         ) %>%
#                         rbind(
#                           prec_recal_T1D_with_T_united_all_genes %>%
#                             mutate(Calculator = "Clinical features and biomarkers")
#                         ),
#                       aes(group = iteration), colour = "grey"
#                     ) + 
#                     geom_path() +
#                     theme_bw() +
#                     facet_grid(~factor(Calculator, levels = c("Clinical features", "Clinical features and biomarkers")), scales = "free",) +
#                     scale_y_continuous("Precision", labels = scales::percent) +
#                     scale_x_continuous("Recall", labels = scales::percent) +
#                     theme_bw() +
#                     geom_label(
#                       data = prec_recal_curves %>%
#                         select(-precision, -recall) %>%
#                         distinct() %>%
#                         mutate(
#                           auc = paste0(" AUC:", signif(auc, 2), " "),
#                           mean = paste0("Mean prob:", signif(mean, 2)*100, "%")
#                         ) %>%
#                         filter(Dataset == "UNITED" & Model == "Type 1") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = c("Biomarkers", "No Biomarkers"), labels = c("Clinical features and biomarkers", "Clinical features"))
#                         ),
#                       mapping = aes(x = -Inf, y = -Inf, label = auc),
#                       size = 7,
#                       label.size = NA,
#                       hjust = -1,
#                       vjust = -6.5
#                     ) +
#                     theme(
#                       panel.spacing.x = unit(1.5, "lines")
#                     )),
#   
#   patchwork::free(prec_recal_curves %>%
#                     filter(Dataset == "UNITED" & Model == "Type 2") %>%
#                     ggplot(aes(x = recall, y = precision)) +
#                     geom_path(
#                       data = prec_recal_T2D_new_united_all_genes %>%
#                         mutate(
#                           Calculator = "No Biomarkers"
#                         ),
#                       aes(group = iteration), colour = "grey"
#                     ) +
#                     geom_path() +
#                     theme_bw() +
#                     scale_y_continuous("Precision", labels = scales::percent) +
#                     scale_x_continuous("Recall", labels = scales::percent) +
#                     theme_bw() +
#                     geom_label(
#                       data = dat_text %>%
#                         filter(Dataset == "UNITED" & Model == "Type 2"),
#                       mapping = aes(x = -Inf, y = -Inf, label = auc),
#                       size = 7,
#                       label.size = NA,
#                       hjust = -0.9,
#                       vjust = -6.5
#                     )),
#   
#   ncol = 1
#   
# ) + patchwork::plot_annotation(tag_levels = list(c("A", "B"))) &
#   theme(
#     axis.text = element_text(size = 14),
#     axis.title = element_text(size = 16),
#     strip.text = element_text(size = 11)
#   )
# 
# pdf("figures/united_prec_recal_thin_100.pdf", width = 8, height = 8)
# plot_prob_prec_recal_united +
#   patchwork::plot_layout(
#     design = "
#     AAAA
#     #BB#
#     "
#   )
# dev.off()
# 
# 
# 
# 
# 
# ################################
# # Testing different roc curves
# 
# 
# plot_prob_rocs_united_testing_1 <- patchwork::wrap_plots(
#   
#   # roc insulin-treated
#   patchwork::free(roc_curves %>%
#                     filter(Dataset == "UNITED" & Model == "Type 1") %>%
#                     mutate(
#                       Calculator = factor(Calculator, levels = c("No Biomarkers", "Biomarkers"), labels = c("Clinical features", "Clinical features and biomarkers")),
#                       iteration = 0
#                     ) %>%
#                     ggplot(aes(x = 1- specificities, y = sensitivities)) +
#                     geom_path() +
#                     theme_bw() +
#                     facet_grid(~factor(Calculator, levels = c("Clinical features", "Clinical features and biomarkers"), labels = c("Early-insulin-treated: clinical features, no biomarkers", "Early-insulin-treated: clinical features, with biomarkers")), scales = "free",labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
#                     scale_y_continuous("Sensitivity", labels = scales::percent) +
#                     scale_x_continuous("1- Specificity", labels = scales::percent) +
#                     theme_bw() +
#                     geom_label(
#                       data = roc_curves %>%
#                         select(-sensitivities, -specificities) %>%
#                         distinct() %>%
#                         mutate(
#                           auc = paste0(" AUC:", signif(auc, 2), " "),
#                           mean = paste0("Mean prob:", signif(mean, 2)*100, "%")
#                         ) %>%
#                         filter(Dataset == "UNITED" & Model == "Type 1") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = c("Biomarkers", "No Biomarkers"), labels = c("Clinical features and biomarkers", "Clinical features"))
#                         ),
#                       mapping = aes(x = -Inf, y = -Inf, label = auc),
#                       size = 7,
#                       label.size = NA,
#                       hjust = -0.8,
#                       vjust = -0.5
#                     ) +
#                     theme(
#                       panel.spacing.x = unit(1.5, "lines")
#                     )),
#   # roc not early insulin treated
#   patchwork::free(roc_curves %>%
#                     filter(Dataset == "UNITED" & Model == "Type 2") %>%
#                     mutate(Calculator = factor(Calculator, levels = "No Biomarkers", labels = "Not-early-insulin-treated")) %>%
#                     ggplot(aes(x = 1- specificities, y = sensitivities)) +
#                     geom_path() +
#                     theme_bw() +
#                     facet_grid(~~factor(Calculator, levels = c("Not-early-insulin-treated"), labels = c("Not-early-insulin-treated: clinical features, no biomarkers")), scales = "free",labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
#                     scale_y_continuous("Sensitivity", labels = scales::percent) +
#                     scale_x_continuous("1- Specificity", labels = scales::percent) +
#                     theme_bw() +
#                     geom_label(
#                       data = dat_text %>%
#                         filter(Dataset == "UNITED" & Model == "Type 2") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = c("Clinical features"), labels = c("Not-early-insulin-treated"))
#                         ),
#                       mapping = aes(x = -Inf, y = -Inf, label = auc),
#                       size = 7,
#                       label.size = NA,
#                       hjust = -0.8,
#                       vjust = -0.5
#                     )),
#   
#   ncol = 1
#   
# ) + patchwork::plot_annotation(tag_levels = list(c("A", "B"))) &
#   theme(
#     axis.text = element_text(size = 14),
#     axis.title = element_text(size = 16),
#     strip.text = element_text(size = 11)
#   )
# 
# #:----------------------
# 
# roc_t1d_no_T <- pROC::roc(response = dataset.UNITED_type1_all_genes$M, predictor = predictions_dataset.UNITED_type1_all_genes_no_T$prob)
# 
# roccoords_t1d_no_T <- pROC::coords(roc_t1d_no_T, ret = c("threshold", "specificity", "sensitivity"), transpose = FALSE)
# 
# roc_t1d_with_T <- pROC::roc(response = dataset.UNITED_type1_all_genes$M, predictor = predictions_dataset.UNITED_type1_all_genes_with_T$prob)
# 
# roccoords_t1d_with_T <- pROC::coords(roc_t1d_with_T, ret = c("threshold", "specificity", "sensitivity"), transpose = FALSE)
# 
# roc_t2d_new <- pROC::roc(response = dataset.UNITED_type2_all_genes$M, predictor = predictions_dataset.UNITED_type2_all_genes_new$prob)
# 
# roccoords_t2d_new <- pROC::coords(roc_t2d_new, ret = c("threshold", "specificity", "sensitivity"), transpose = FALSE)
# 
# 
# roc_curves_new <- data.frame(
#   ci.coords(roc_t1d_no_T, roccoords_t1d_no_T$threshold, ret= c("specificity", "sensitivity"))
# ) %>%
#   mutate(
#     Dataset = "UNITED", Model = "Type 1", Calculator = "No Biomarkers"
#   ) %>%
#   rbind(
#     data.frame(
#       ci.coords(roc_t1d_with_T, roccoords_t1d_with_T$threshold, ret= c("specificity", "sensitivity"))
#     ) %>%
#       mutate(
#         Dataset = "UNITED", Model = "Type 1", Calculator = "Biomarkers"
#       ),
#     data.frame(
#       ci.coords(roc_t2d_new, roccoords_t2d_new$threshold, ret= c("specificity", "sensitivity"))
#     ) %>%
#       mutate(
#         Dataset = "UNITED", Model = "Type 2", Calculator = "No Biomarkers"
#       )
#   )
# 
# 
# plot_prob_rocs_united_testing_2 <- patchwork::wrap_plots(
#   
#   # roc insulin-treated
#   patchwork::free(roc_curves_new %>%
#                     filter(Dataset == "UNITED" & Model == "Type 1") %>%
#                     mutate(
#                       Calculator = factor(Calculator, levels = c("No Biomarkers", "Biomarkers"), labels = c("Clinical features", "Clinical features and biomarkers")),
#                       iteration = 0
#                     ) %>%
#                     ggplot() +
#                     geom_path(aes(x = 1- `specificity.2.5.`, y = `sensitivity.2.5.`), linetype = "dashed", colour = "grey") +
#                     geom_path(aes(x = 1- `specificity.97.5.`, y = `sensitivity.97.5.`), linetype = "dashed", colour = "grey") +
#                     geom_path(aes(x = 1- `specificity.50.`, y = `sensitivity.50.`), colour = "black") +
#                     theme_bw() +
#                     facet_grid(~factor(Calculator, levels = c("Clinical features", "Clinical features and biomarkers"), labels = c("Early-insulin-treated: clinical features, no biomarkers", "Early-insulin-treated: clinical features, with biomarkers")), scales = "free",labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
#                     scale_y_continuous("Sensitivity", labels = scales::percent) +
#                     scale_x_continuous("1- Specificity", labels = scales::percent) +
#                     theme_bw() +
#                     geom_label(
#                       data = roc_curves %>%
#                         select(-sensitivities, -specificities) %>%
#                         distinct() %>%
#                         mutate(
#                           auc = paste0(" AUC:", signif(auc, 2), " "),
#                           mean = paste0("Mean prob:", signif(mean, 2)*100, "%")
#                         ) %>%
#                         filter(Dataset == "UNITED" & Model == "Type 1") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = c("Biomarkers", "No Biomarkers"), labels = c("Clinical features and biomarkers", "Clinical features"))
#                         ),
#                       mapping = aes(x = -Inf, y = -Inf, label = auc),
#                       size = 7,
#                       label.size = NA,
#                       hjust = -0.8,
#                       vjust = -0.5
#                     ) +
#                     theme(
#                       panel.spacing.x = unit(1.5, "lines")
#                     )),
#   # roc not early insulin treated
#   patchwork::free(roc_curves_new %>%
#                     filter(Dataset == "UNITED" & Model == "Type 2") %>%
#                     mutate(Calculator = factor(Calculator, levels = "No Biomarkers", labels = "Not-early-insulin-treated")) %>%
#                     ggplot() +
#                     geom_path(aes(x = 1- `specificity.2.5.`, y = `sensitivity.2.5.`), linetype = "dashed", colour = "grey") +
#                     geom_path(aes(x = 1- `specificity.97.5.`, y = `sensitivity.97.5.`), linetype = "dashed", colour = "grey") +
#                     geom_path(aes(x = 1- `specificity.50.`, y = `sensitivity.50.`), colour = "black") +
#                     theme_bw() +
#                     facet_grid(~factor(Calculator, levels = c("Not-early-insulin-treated"), labels = c("Not-early-insulin-treated: clinical features, no biomarkers")), scales = "free",labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
#                     scale_y_continuous("Sensitivity", labels = scales::percent) +
#                     scale_x_continuous("1- Specificity", labels = scales::percent) +
#                     theme_bw() +
#                     geom_label(
#                       data = dat_text %>%
#                         filter(Dataset == "UNITED" & Model == "Type 2") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = c("Clinical features"), labels = c("Not-early-insulin-treated"))
#                         ),
#                       mapping = aes(x = -Inf, y = -Inf, label = auc),
#                       size = 7,
#                       label.size = NA,
#                       hjust = -0.8,
#                       vjust = -0.5
#                     )),
#   
#   ncol = 1
#   
# ) + patchwork::plot_annotation(tag_levels = list(c("A", "B"))) &
#   theme(
#     axis.text = element_text(size = 14),
#     axis.title = element_text(size = 16),
#     strip.text = element_text(size = 11)
#   )
# 
# 
# #:----------------------
# 
# 
# roc_curves_new_polygon <- data.frame(
#   x = c(
#     ci.coords(roc_t1d_no_T, roccoords_t1d_no_T$threshold, ret= c("specificity", "sensitivity"))$specificity[,1], 
#     rev(ci.coords(roc_t1d_no_T, roccoords_t1d_no_T$threshold, ret= c("specificity", "sensitivity"))$specificity[,3])
#   ),
#   y = c(
#     ci.coords(roc_t1d_no_T, roccoords_t1d_no_T$threshold, ret= c("specificity", "sensitivity"))$sensitivity[,1], 
#     rev(ci.coords(roc_t1d_no_T, roccoords_t1d_no_T$threshold, ret= c("specificity", "sensitivity"))$sensitivity[,3])
#   )
# ) %>%
#   mutate(
#     Dataset = "UNITED", Model = "Type 1", Calculator = "No Biomarkers"
#   ) %>%
#   rbind(
#     data.frame(
#       x = c(
#         ci.coords(roc_t1d_with_T, roccoords_t1d_with_T$threshold, ret= c("specificity", "sensitivity"))$specificity[,1], 
#         rev(ci.coords(roc_t1d_with_T, roccoords_t1d_with_T$threshold, ret= c("specificity", "sensitivity"))$specificity[,3])
#       ),
#       y = c(
#         ci.coords(roc_t1d_with_T, roccoords_t1d_with_T$threshold, ret= c("specificity", "sensitivity"))$sensitivity[,1], 
#         rev(ci.coords(roc_t1d_with_T, roccoords_t1d_with_T$threshold, ret= c("specificity", "sensitivity"))$sensitivity[,3])
#       )
#     ) %>%
#       mutate(
#         Dataset = "UNITED", Model = "Type 1", Calculator = "Biomarkers"
#       ),
#     data.frame(
#       x = c(
#         ci.coords(roc_t2d_new, roccoords_t2d_new$threshold, ret= c("specificity", "sensitivity"))$specificity[,1], 
#         rev(ci.coords(roc_t2d_new, roccoords_t2d_new$threshold, ret= c("specificity", "sensitivity"))$specificity[,3])
#       ),
#       y = c(
#         ci.coords(roc_t2d_new, roccoords_t2d_new$threshold, ret= c("specificity", "sensitivity"))$sensitivity[,1], 
#         rev(ci.coords(roc_t2d_new, roccoords_t2d_new$threshold, ret= c("specificity", "sensitivity"))$sensitivity[,3])
#       )
#     ) %>%
#       mutate(
#         Dataset = "UNITED", Model = "Type 2", Calculator = "No Biomarkers"
#       )
#   )
# 
# 
# 
# plot_prob_rocs_united_testing_3 <- patchwork::wrap_plots(
#   
#   # roc insulin-treated
#   patchwork::free(roc_curves_new %>%
#                     filter(Dataset == "UNITED" & Model == "Type 1") %>%
#                     mutate(
#                       Calculator = factor(Calculator, levels = c("No Biomarkers", "Biomarkers"), labels = c("Clinical features", "Clinical features and biomarkers"))
#                     ) %>%
#                     ggplot() +
#                     geom_polygon(
#                       data = roc_curves_new_polygon %>%
#                         filter(Dataset == "UNITED" & Model == "Type 1") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = c("No Biomarkers", "Biomarkers"), labels = c("Clinical features", "Clinical features and biomarkers"))
#                         ),
#                       aes(x = 1 - x, y = y), fill = "grey") +
#                     geom_path(aes(x = 1- `specificity.50.`, y = `sensitivity.50.`), colour = "black") +
#                     theme_bw() +
#                     facet_grid(~factor(Calculator, levels = c("Clinical features", "Clinical features and biomarkers"), labels = c("Early-insulin-treated: clinical features, no biomarkers", "Early-insulin-treated: clinical features, with biomarkers")), scales = "free",labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
#                     scale_y_continuous("Sensitivity", labels = scales::percent) +
#                     scale_x_continuous("1- Specificity", labels = scales::percent) +
#                     theme_bw() +
#                     geom_label(
#                       data = roc_curves %>%
#                         select(-sensitivities, -specificities) %>%
#                         distinct() %>%
#                         mutate(
#                           auc = paste0(" AUC:", signif(auc, 2), " "),
#                           mean = paste0("Mean prob:", signif(mean, 2)*100, "%")
#                         ) %>%
#                         filter(Dataset == "UNITED" & Model == "Type 1") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = c("Biomarkers", "No Biomarkers"), labels = c("Clinical features and biomarkers", "Clinical features"))
#                         ),
#                       mapping = aes(x = -Inf, y = -Inf, label = auc),
#                       size = 7,
#                       label.size = NA,
#                       hjust = -0.8,
#                       vjust = -0.5
#                     ) +
#                     theme(
#                       panel.spacing.x = unit(1.5, "lines")
#                     )),
#   # roc not early insulin treated
#   patchwork::free(roc_curves_new %>%
#                     filter(Dataset == "UNITED" & Model == "Type 2") %>%
#                     mutate(Calculator = factor(Calculator, levels = "No Biomarkers", labels = "Not-early-insulin-treated")) %>%
#                     ggplot() +
#                     geom_polygon(
#                       data = roc_curves_new_polygon %>%
#                         filter(Dataset == "UNITED" & Model == "Type 2") %>%
#                         mutate(Calculator = factor(Calculator, levels = "No Biomarkers", labels = "Not-early-insulin-treated")),
#                       aes(x = 1 - x, y = y), fill = "grey") +
#                     geom_path(aes(x = 1- `specificity.50.`, y = `sensitivity.50.`), colour = "black") +
#                     theme_bw() +
#                     facet_grid(~factor(Calculator, levels = c("Not-early-insulin-treated"), labels = c("Not-early-insulin-treated: clinical features, no biomarkers")), scales = "free",labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
#                     scale_y_continuous("Sensitivity", labels = scales::percent) +
#                     scale_x_continuous("1- Specificity", labels = scales::percent) +
#                     theme_bw() +
#                     geom_label(
#                       data = dat_text %>%
#                         filter(Dataset == "UNITED" & Model == "Type 2") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = c("Clinical features"), labels = c("Not-early-insulin-treated"))
#                         ),
#                       mapping = aes(x = -Inf, y = -Inf, label = auc),
#                       size = 7,
#                       label.size = NA,
#                       hjust = -0.8,
#                       vjust = -0.5
#                     )),
#   
#   ncol = 1
#   
# ) + patchwork::plot_annotation(tag_levels = list(c("A", "B"))) &
#   theme(
#     axis.text = element_text(size = 14),
#     axis.title = element_text(size = 16),
#     strip.text = element_text(size = 11)
#   )
# 
# 
# 
# #:----------------------
# 
# 
# 
# ### No biomarker models
# auc_T1D_no_T_united_all_genes <- calc_auroc(dataset.UNITED_type1_all_genes$M, predictions_dataset.UNITED_type1_all_genes_no_T_full, thinning = 10)
# 
# min_no_T <- order(auc_T1D_no_T_united_all_genes)[0.025*length(auc_T1D_no_T_united_all_genes)]
# max_no_T <- order(auc_T1D_no_T_united_all_genes)[0.975*length(auc_T1D_no_T_united_all_genes)]
# 
# coords_min_no_T <- pROC::coords(pROC::roc(response = dataset.UNITED_type1_all_genes$M, predictor = predictions_dataset.UNITED_type1_all_genes_no_T_full[seq(1, 800000, 10)[min_no_T],]))
# coords_max_no_T <- pROC::coords(pROC::roc(response = dataset.UNITED_type1_all_genes$M, predictor = predictions_dataset.UNITED_type1_all_genes_no_T_full[seq(1, 800000, 10)[max_no_T],]))
# 
# d1 <- coords_min_no_T %>%
#   mutate(type = "min") %>%
#   rbind(
#     coords_max_no_T %>%
#       mutate(type = "max")
#   ) %>%
#   mutate(
#     Dataset = "UNITED", Model = "Type 1", Calculator = "No Biomarkers"
#   )
# 
# ### Biomarker models
# auc_T1D_with_T_united_all_genes <- calc_auroc(dataset.UNITED_type1_all_genes$M, predictions_dataset.UNITED_type1_all_genes_with_T_full, thinning = 10)
# 
# min_with_T <- order(auc_T1D_with_T_united_all_genes)[0.025*length(auc_T1D_with_T_united_all_genes)]
# max_with_T <- order(auc_T1D_with_T_united_all_genes)[0.975*length(auc_T1D_with_T_united_all_genes)]
# 
# coords_min_with_T <- pROC::coords(pROC::roc(response = dataset.UNITED_type1_all_genes$M, predictor = predictions_dataset.UNITED_type1_all_genes_with_T_full[seq(1, 800000, 10)[min_with_T],]))
# coords_max_with_T <- pROC::coords(pROC::roc(response = dataset.UNITED_type1_all_genes$M, predictor = predictions_dataset.UNITED_type1_all_genes_with_T_full[seq(1, 800000, 10)[max_with_T],]))
# 
# d2 <- coords_min_with_T %>%
#   mutate(type = "min") %>%
#   rbind(
#     coords_max_with_T %>%
#       mutate(type = "max")
#   ) %>%
#   mutate(
#     Dataset = "UNITED", Model = "Type 1", Calculator = "Biomarkers"
#   )
# 
# ## Type 2 UNITED
# auc_T2D_new_united_all_genes <- calc_auroc(dataset.UNITED_type2_all_genes$M, predictions_dataset.UNITED_type2_all_genes_new_full, thinning = 10)
# 
# min_new <- order(auc_T2D_new_united_all_genes)[0.025*length(auc_T2D_new_united_all_genes)]
# max_new <- order(auc_T2D_new_united_all_genes)[0.975*length(auc_T2D_new_united_all_genes)]
# 
# coords_min_new <- pROC::coords(pROC::roc(response = dataset.UNITED_type2_all_genes$M, predictor = predictions_dataset.UNITED_type2_all_genes_new_full[seq(1, 800000, 10)[min_new],]))
# coords_max_new <- pROC::coords(pROC::roc(response = dataset.UNITED_type2_all_genes$M, predictor = predictions_dataset.UNITED_type2_all_genes_new_full[seq(1, 800000, 10)[max_new],]))
# 
# d3 <- coords_min_new %>%
#   mutate(type = "min") %>%
#   rbind(
#     coords_max_new %>%
#       mutate(type = "max")
#   ) %>%
#   mutate(
#     Dataset = "UNITED", Model = "Type 2", Calculator = "No Biomarkers"
#   )
# 
# 
# 
# roc_dataset_test <- rbind(d1, d2, d3)
# 
# 
# plot_prob_rocs_united_testing_4 <- patchwork::wrap_plots(
#   
#   # roc insulin-treated
#   patchwork::free(roc_curves %>%
#                     filter(Dataset == "UNITED" & Model == "Type 1") %>%
#                     mutate(
#                       Calculator = factor(Calculator, levels = c("No Biomarkers", "Biomarkers"), labels = c("Clinical features", "Clinical features and biomarkers"))
#                     ) %>%
#                     ggplot() +
#                     geom_path(
#                       data = roc_dataset_test %>%
#                         filter(type == "min") %>%
#                         filter(Dataset == "UNITED" & Model == "Type 1") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = c("No Biomarkers", "Biomarkers"), labels = c("Clinical features", "Clinical features and biomarkers"))
#                         ),
#                       aes(x = 1 - `specificity`, y = `sensitivity`), colour = "blue"
#                     ) +
#                     geom_path(
#                       data = roc_dataset_test %>%
#                         filter(type == "max") %>%
#                         filter(Dataset == "UNITED" & Model == "Type 1") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = c("No Biomarkers", "Biomarkers"), labels = c("Clinical features", "Clinical features and biomarkers"))
#                         ),
#                       aes(x = 1 - `specificity`, y = `sensitivity`), colour = "red"
#                     ) +
#                     geom_path(aes(x = 1- specificities, y = sensitivities)) +
#                     theme_bw() +
#                     facet_grid(~factor(Calculator, levels = c("Clinical features", "Clinical features and biomarkers"), labels = c("Early-insulin-treated: clinical features, no biomarkers", "Early-insulin-treated: clinical features, with biomarkers")), scales = "free",labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
#                     scale_y_continuous("Sensitivity", labels = scales::percent) +
#                     scale_x_continuous("1- Specificity", labels = scales::percent) +
#                     theme_bw() +
#                     geom_label(
#                       data = roc_curves %>%
#                         select(-sensitivities, -specificities) %>%
#                         distinct() %>%
#                         mutate(
#                           auc = paste0(" AUC:", signif(auc, 2), " "),
#                           mean = paste0("Mean prob:", signif(mean, 2)*100, "%")
#                         ) %>%
#                         filter(Dataset == "UNITED" & Model == "Type 1") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = c("Biomarkers", "No Biomarkers"), labels = c("Clinical features and biomarkers", "Clinical features"))
#                         ),
#                       mapping = aes(x = -Inf, y = -Inf, label = auc),
#                       size = 7,
#                       label.size = NA,
#                       hjust = -0.8,
#                       vjust = -0.5
#                     ) +
#                     theme(
#                       panel.spacing.x = unit(1.5, "lines")
#                     )),
#   # roc not early insulin treated
#   patchwork::free(roc_curves %>%
#                     filter(Dataset == "UNITED" & Model == "Type 2") %>%
#                     mutate(Calculator = factor(Calculator, levels = "No Biomarkers", labels = "Not-early-insulin-treated")) %>%
#                     ggplot() +
#                     geom_path(
#                       data = roc_dataset_test %>%
#                         filter(type == "min") %>%
#                         filter(Dataset == "UNITED" & Model == "Type 2") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = "No Biomarkers", labels = "Not-early-insulin-treated")
#                         ),
#                       aes(x = 1 - `specificity`, y = `sensitivity`), colour = "blue"
#                     ) +
#                     geom_path(
#                       data = roc_dataset_test %>%
#                         filter(type == "max") %>%
#                         filter(Dataset == "UNITED" & Model == "Type 2") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = "No Biomarkers", labels = "Not-early-insulin-treated")
#                         ),
#                       aes(x = 1 - `specificity`, y = `sensitivity`), colour = "red"
#                     ) +
#                     geom_path(aes(x = 1- specificities, y = sensitivities)) +
#                     theme_bw() +
#                     facet_grid(~factor(Calculator, levels = c("Not-early-insulin-treated"), labels = c("Not-early-insulin-treated: clinical features, no biomarkers")), scales = "free",labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
#                     scale_y_continuous("Sensitivity", labels = scales::percent) +
#                     scale_x_continuous("1- Specificity", labels = scales::percent) +
#                     theme_bw() +
#                     geom_label(
#                       data = dat_text %>%
#                         filter(Dataset == "UNITED" & Model == "Type 2") %>%
#                         mutate(
#                           Calculator = factor(Calculator, levels = c("Clinical features"), labels = c("Not-early-insulin-treated"))
#                         ),
#                       mapping = aes(x = -Inf, y = -Inf, label = auc),
#                       size = 7,
#                       label.size = NA,
#                       hjust = -0.8,
#                       vjust = -0.5
#                     )),
#   
#   ncol = 1
#   
# ) + patchwork::plot_annotation(tag_levels = list(c("A", "B"))) &
#   theme(
#     axis.text = element_text(size = 14),
#     axis.title = element_text(size = 16),
#     strip.text = element_text(size = 11)
#   )
# 
# 
# 
# 
# 
# 
# pdf("figures/united_roc_thin_100_testing.pdf", width = 9, height = 8)
# plot_prob_rocs_united +
#   patchwork::plot_layout(
#     design = "
#     AAAA
#     #BB#
#     "
#   )
# plot_prob_rocs_united_testing_1 +
#   patchwork::plot_layout(
#     design = "
#     AAAA
#     #BB#
#     "
#   )
# plot_prob_rocs_united_testing_2 +
#   patchwork::plot_layout(
#     design = "
#     AAAA
#     #BB#
#     "
#   )
# plot_prob_rocs_united_testing_3 +
#   patchwork::plot_layout(
#     design = "
#     AAAA
#     #BB#
#     "
#   )
# plot_prob_rocs_united_testing_4 +
#   patchwork::plot_layout(
#     design = "
#     AAAA
#     #BB#
#     "
#   )
# 
# dev.off()














