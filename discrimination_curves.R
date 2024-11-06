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
predictions_dataset.UNITED_type1_all_genes_no_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_all_genes_no_T_full.rds")
predictions_dataset.UNITED_type1_all_genes_with_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_all_genes_with_T_full.rds")
predictions_dataset.UNITED_type2_all_genes_new <- readRDS("model_predictions/predictions_dataset.UNITED_type2_all_genes_new_full.rds")

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
auc_T1D_no_T_united_all_genes <- calc_auroc(dataset.UNITED_type1_all_genes$M, predictions_dataset.UNITED_type1_all_genes_no_T, thinning = 10)
# quantile(auc_T1D_no_T_united_all_genes, probs = c(0.025, 0.5, 0.975)) # thinning = 10
# # 2.5%       50%     97.5% 
# # 0.7986016 0.8357906 0.8562092 

### Biomarker models
auc_T1D_with_T_united_all_genes <- calc_auroc(dataset.UNITED_type1_all_genes$M, predictions_dataset.UNITED_type1_all_genes_with_T, thinning = 10)
# quantile(auc_T1D_with_T_united_all_genes, probs = c(0.025, 0.5, 0.975)) # thinning = 10
# # 2.5%       50%     97.5% 
# # 0.9765922 0.9792268 0.9812028 

## Type 2 UNITED
auc_T2D_new_united_all_genes <- calc_auroc(dataset.UNITED_type2_all_genes$M, predictions_dataset.UNITED_type2_all_genes_new, thinning = 10)
# quantile(auc_T2D_new_united_all_genes, probs = c(0.025, 0.5, 0.975)) # thinning = 10
# # 2.5%       50%     97.5% 
# # 0.7892308 0.8042308 0.8192308 


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
pr_auc_T1D_no_T_united_all_genes <- calc_auc_pr(dataset.UNITED_type1_all_genes$M, predictions_dataset.UNITED_type1_all_genes_no_T, thinning = 10)
# quantile(pr_auc_T1D_no_T_united_all_genes, probs = c(0.025, 0.5, 0.975)) # thinning = 10
# 2.5%        50%      97.5% 
# 0.02233093 0.03205236 0.04040395 

### Biomarker models
pr_auc_T1D_with_T_united_all_genes <- calc_auc_pr(dataset.UNITED_type1_all_genes$M, predictions_dataset.UNITED_type1_all_genes_with_T, thinning = 10)
# quantile(pr_auc_T1D_with_T_united_all_genes, probs = c(0.025, 0.5, 0.975)) # thinning = 10
# 2.5%       50%     97.5% 
# 0.1705823 0.2104027 0.2396283 

## Type 2 UNITED
pr_auc_T2D_new_united_all_genes <- calc_auc_pr(dataset.UNITED_type2_all_genes$M, predictions_dataset.UNITED_type2_all_genes_new, thinning = 10)
# quantile(pr_auc_T2D_new_united_all_genes, probs = c(0.025, 0.5, 0.975)) # thinning = 10
# # 2.5%       50%     97.5% 
# # 0.5601752 0.5872483 0.6154755 

#:------------------------------------------------------------

# Calculate ROC with intervals
calc_roc <- function(data, predictions, thinning = 100) {
  
  output <- NULL
  
  # sequence to iterate through
  sequence_list <- seq(1, nrow(predictions), thinning)
  
  for (i in 1:length(sequence_list)) {
    
    # print the current iteration
    if (i %% 1000 == 0) {
      print(paste(i, "out of", length(sequence_list)))
    }
    
    ## calculate ROC
    interim <- pROC::roc(response = data, predictor = predictions[sequence_list[i],], levels = c(0,1), direction = "<") %>%
      magrittr::extract(c(2:3)) %>%
      as.data.frame() %>%
      mutate(iteration = paste0(sequence_list[i]))
    
    output <- rbind(output, interim)
    
  }
  
  return(output)
  
}


## Type 1 UNITED

### No biomarker models
roc_T1D_no_T_united_all_genes <- calc_roc(dataset.UNITED_type1_all_genes$M, predictions_dataset.UNITED_type1_all_genes_no_T, thinning = 100)

# plot for ROC with grey being iterations, black being the ROC for average prediction
plot_roc_T1D_no_T_united_all_genes <- ggplot() +
  ## all iterations
  geom_path(
    data = roc_T1D_no_T_united_all_genes,
    aes(x = 1-specificities, y= sensitivities, group = iteration), colour = "grey"
  ) +
  ## average predictions
  geom_path(
    data = pROC::roc(response = dataset.UNITED_type1_all_genes$M, predictor = colMeans(predictions_dataset.UNITED_type1_all_genes_no_T), levels = c(0,1), direction = "<") %>%
      magrittr::extract(c(2:3)) %>%
      as.data.frame(),
    aes(x = 1-specificities, y= sensitivities), colour = "black"
  )

### Biomarker models
roc_T1D_with_T_united_all_genes <- calc_roc(dataset.UNITED_type1_all_genes$M, predictions_dataset.UNITED_type1_all_genes_with_T, thinning = 100)

# plot for ROC with grey being iterations, black being the ROC for average prediction
plot_roc_T1D_with_T_united_all_genes <- ggplot() +
  ## all iterations
  geom_path(
    data = roc_T1D_with_T_united_all_genes,
    aes(x = 1-specificities, y= sensitivities, group = iteration), colour = "grey"
  ) +
  ## average predictions
  geom_path(
    data = pROC::roc(response = dataset.UNITED_type1_all_genes$M, predictor = colMeans(predictions_dataset.UNITED_type1_all_genes_with_T)) %>%
      magrittr::extract(c(2:3)) %>%
      as.data.frame(),
    aes(x = 1-specificities, y= sensitivities), colour = "black"
  )

## Type 2 UNITED
roc_T2D_new_united_all_genes <- calc_roc(dataset.UNITED_type2_all_genes$M, predictions_dataset.UNITED_type2_all_genes_new, thinning = 100)

# plot for ROC with grey being iterations, black being the ROC for average prediction
plot_roc_T2D_new_united_all_genes <- ggplot() +
  ## all iterations
  geom_path(
    data = roc_T2D_new_united_all_genes,
    aes(x = 1-specificities, y= sensitivities, group = iteration), colour = "grey"
  ) +
  ## average predictions
  geom_path(
    data = pROC::roc(response = dataset.UNITED_type2_all_genes$M, predictor = colMeans(predictions_dataset.UNITED_type2_all_genes_new)) %>%
      magrittr::extract(c(2:3)) %>%
      as.data.frame(),
    aes(x = 1-specificities, y= sensitivities), colour = "black"
  )



#:------------------------------------------------------------


# Calculate Precision-recall with intervals
calc_prec_recal_curve <- function(data, predictions, thinning = 100) {
  
  output <- NULL
  
  # sequence to iterate through
  sequence_list <- seq(1, nrow(predictions), thinning)
  
  for (i in 1:length(sequence_list)) {
    
    # print the current iteration
    if (i %% 100 == 0) {
      print(paste(i, "out of", length(sequence_list)))
    }
    
    ## calculate ROC
    interim <- pROC::coords(pROC::roc(response = data, predictor = predictions[sequence_list[i],], levels = c(0,1), direction = "<"), ret = c("precision", "recall")) %>%
      as.data.frame() %>%
      mutate(iteration = paste0(sequence_list[i]))
    
    output <- rbind(output, interim)
    
  }
  
  return(output)
  
}


## Type 1 UNITED

### No biomarker models
prec_recal_T1D_no_T_united_all_genes <- calc_prec_recal_curve(dataset.UNITED_type1_all_genes$M, predictions_dataset.UNITED_type1_all_genes_no_T, thinning = 100)

# plot for ROC with grey being iterations, black being the ROC for average prediction
plot_prec_recal_T1D_no_T_united_all_genes <- ggplot() +
  ## all iterations
  geom_path(
    data = prec_recal_T1D_no_T_united_all_genes,
    aes(x = recall, y= precision, group = iteration), colour = "grey"
  ) +
  ## average predictions
  geom_path(
    data = pROC::coords(pROC::roc(response = dataset.UNITED_type1_all_genes$M, predictor = colMeans(predictions_dataset.UNITED_type1_all_genes_no_T)), ret = c("precision", "recall")) %>%
      as.data.frame(),
    aes(x = recall, y= precision), colour = "black"
  )

### Biomarker models
prec_recal_T1D_with_T_united_all_genes <- calc_prec_recal_curve(dataset.UNITED_type1_all_genes$M, predictions_dataset.UNITED_type1_all_genes_with_T, thinning = 100)

# plot for ROC with grey being iterations, black being the ROC for average prediction
plot_prec_recal_T1D_with_T_united_all_genes <- ggplot() +
  ## all iterations
  geom_path(
    data = prec_recal_T1D_with_T_united_all_genes,
    aes(x = recall, y= precision, group = iteration), colour = "grey"
  ) +
  ## average predictions
  geom_path(
    data = pROC::coords(pROC::roc(response = dataset.UNITED_type1_all_genes$M, predictor = colMeans(predictions_dataset.UNITED_type1_all_genes_with_T)), ret = c("precision", "recall")) %>%
      as.data.frame(),
    aes(x = recall, y= precision), colour = "black"
  )

## Type 2 UNITED
prec_recal_T2D_new_united_all_genes <- calc_prec_recal_curve(dataset.UNITED_type2_all_genes$M, predictions_dataset.UNITED_type2_all_genes_new, thinning = 100)

# plot for ROC with grey being iterations, black being the ROC for average prediction
plot_prec_recal_T2D_new_united_all_genes <- ggplot() +
  ## all iterations
  geom_path(
    data = prec_recal_T2D_new_united_all_genes,
    aes(x = recall, y= precision, group = iteration), colour = "grey"
  ) +
  ## average predictions
  geom_path(
    data = pROC::coords(pROC::roc(response = dataset.UNITED_type2_all_genes$M, predictor = colMeans(predictions_dataset.UNITED_type2_all_genes_new)), ret = c("precision", "recall")) %>%
      as.data.frame(),
    aes(x = recall, y= precision), colour = "black"
  )


#:--------------------------------------------------

# Boxplot and roc curves

roc_curves <- data.frame(prob = colMeans(predictions_dataset.UNITED_type1_all_genes_with_T)) %>%
  cbind(Mody = dataset.UNITED_type1_all_genes$M) %>%
  pROC::roc(response = Mody, predictor = prob) %>%
  magrittr::extract(2:3) %>%
  as.data.frame() %>%
  mutate(
    ROCAUC =  unname(data.frame(prob = colMeans(predictions_dataset.UNITED_type1_all_genes_with_T)) %>%
                       cbind(Mody = dataset.UNITED_type1_all_genes$M) %>%
                       pROC::roc(response = Mody, predictor = prob) %>%
                       magrittr::extract(c(9)) %>%
                       unlist()),
    mean = mean(colMeans(predictions_dataset.UNITED_type1_all_genes_with_T), na.rm = TRUE)
  ) %>%
  mutate(Dataset = "UNITED", Model = "Type 1", Calculator = "Biomarkers") %>%
  rbind(
    data.frame(prob = colMeans(predictions_dataset.UNITED_type1_all_genes_no_T)) %>%
      cbind(Mody = dataset.UNITED_type1_all_genes$M) %>%
      pROC::roc(response = Mody, predictor = prob) %>%
      magrittr::extract(2:3) %>%
      as.data.frame() %>%
      mutate(
        ROCAUC = unname(data.frame(prob = colMeans(predictions_dataset.UNITED_type1_all_genes_no_T)) %>%
                          cbind(Mody = dataset.UNITED_type1_all_genes$M) %>%
                          pROC::roc(response = Mody, predictor = prob) %>%
                          magrittr::extract(c(9)) %>%
                          unlist()),
        mean = mean(colMeans(predictions_dataset.UNITED_type1_all_genes_no_T), na.rm = TRUE)
      ) %>%
      mutate(Dataset = "UNITED", Model = "Type 1", Calculator = "No Biomarkers"), 
    data.frame(prob = colMeans(predictions_dataset.UNITED_type2_all_genes_new)) %>%
      cbind(Mody = dataset.UNITED_type2_all_genes$M) %>%
      pROC::roc(response = Mody, predictor = prob) %>%
      magrittr::extract(2:3) %>%
      as.data.frame() %>%
      mutate(
        ROCAUC = unname(data.frame(prob = colMeans(predictions_dataset.UNITED_type2_all_genes_new)) %>%
                          cbind(Mody = dataset.UNITED_type2_all_genes$M) %>%
                          pROC::roc(response = Mody, predictor = prob) %>%
                          magrittr::extract(c(9)) %>%
                          unlist()),
        mean = mean(colMeans(predictions_dataset.UNITED_type2_all_genes_new), na.rm = TRUE)
      ) %>%
      mutate(Dataset = "UNITED", Model = "Type 2", Calculator = " ")
  ) %>%
  rename("auc" = "ROCAUC")


dat_text <- roc_curves %>%
  select(-sensitivities, -specificities) %>%
  distinct() %>%
  mutate(
    auc = paste0(" AUC:", signif(auc, 2), " "),
    mean = paste0("Mean prob:", signif(mean, 2)*100, "%"),
    Calculator = factor(Calculator, levels = c("Biomarkers", "No Biomarkers"), labels = c("Clinical features and biomarkers", "Clinical features"))
  )



plot_prob_boxplot_rocs_united <- patchwork::wrap_plots(
  
  # Panel A - UNITED insulin-treated
  
  patchwork::wrap_plots(
    
    # Boxplots
    dataset.UNITED_type1_all_genes %>%
      select(M) %>%
      rename("Mody" = "M") %>%
      cbind(
        prob_with = colMeans(predictions_dataset.UNITED_type1_all_genes_with_T),
        prob_without = colMeans(predictions_dataset.UNITED_type1_all_genes_no_T)
      ) %>%
      gather("key", "Probability", -Mody) %>% 
      mutate(
        Mody = factor(Mody, levels = c(0, 1), labels = c("Negative", "Positive")),
        key = factor(key, levels = c("prob_without", "prob_with"), labels = c("Clinical features", "Clinical features and biomarkers"))
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
        Calculator = factor(Calculator, levels = c("No Biomarkers", "Biomarkers"), labels = c("Clinical features", "Clinical features and biomarkers")),
        iteration = 0
      ) %>%
      ggplot(aes(x = 1- specificities, y = sensitivities)) +
      geom_path(
        data = roc_T1D_no_T_united_all_genes %>%
          mutate(
            Calculator = "Clinical features"
          ) %>%
          rbind(
            roc_T1D_with_T_united_all_genes %>%
              mutate(Calculator = "Clinical features and biomarkers")
          ),
        aes(group = iteration), colour = "grey"
      ) + 
      geom_path() +
      theme_bw() +
      facet_grid(~factor(Calculator, levels = c("Clinical features", "Clinical features and biomarkers")), scales = "free",) +
      scale_y_continuous("Sensitivity", labels = scales::percent) +
      scale_x_continuous("1- Specificity", labels = scales::percent) +
      theme_bw() +
      geom_label(
        data = roc_curves %>%
          select(-sensitivities, -specificities) %>%
          distinct() %>%
          mutate(
            auc = paste0(" AUC:", signif(auc, 2), " "),
            mean = paste0("Mean prob:", signif(mean, 2)*100, "%")
          ) %>%
          filter(Dataset == "UNITED" & Model == "Type 1") %>%
          mutate(
            Calculator = factor(Calculator, levels = c("Biomarkers", "No Biomarkers"), labels = c("Clinical features and biomarkers", "Clinical features"))
          ),
        mapping = aes(x = -Inf, y = -Inf, label = auc),
        size = 7,
        label.size = NA,
        hjust = -0.4,
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
    dataset.UNITED_type2_all_genes %>%
      select(M) %>%
      rename("Mody" = "M") %>%
      cbind(
        Probability = colMeans(predictions_dataset.UNITED_type2_all_genes_new),
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
        data = roc_T2D_new_united_all_genes %>%
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
        mapping = aes(x = -Inf, y = -Inf, label = auc),
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
    strip.text = element_text(size = 11)
  )


#:-------------------------------------------------------------
# Making plots
pdf("figures/united_boxplot_roc_thin_100.pdf", width = 13, height = 9)
plot_prob_boxplot_rocs_united
dev.off()



plot_prob_rocs_united <- patchwork::wrap_plots(
  
  # roc insulin-treated
  free(roc_curves %>%
    filter(Dataset == "UNITED" & Model == "Type 1") %>%
    mutate(
      Calculator = factor(Calculator, levels = c("No Biomarkers", "Biomarkers"), labels = c("Clinical features", "Clinical features and biomarkers")),
      iteration = 0
    ) %>%
    ggplot(aes(x = 1- specificities, y = sensitivities)) +
    geom_path(
      data = roc_T1D_no_T_united_all_genes %>%
        mutate(
          Calculator = "Clinical features"
        ) %>%
        rbind(
          roc_T1D_with_T_united_all_genes %>%
            mutate(Calculator = "Clinical features and biomarkers")
        ),
      aes(group = iteration), colour = "grey"
    ) + 
    geom_path() +
    theme_bw() +
    facet_grid(~factor(Calculator, levels = c("Clinical features", "Clinical features and biomarkers")), scales = "free",) +
    scale_y_continuous("Sensitivity", labels = scales::percent) +
    scale_x_continuous("1- Specificity", labels = scales::percent) +
    theme_bw() +
    geom_label(
      data = roc_curves %>%
        select(-sensitivities, -specificities) %>%
        distinct() %>%
        mutate(
          auc = paste0(" AUC:", signif(auc, 2), " "),
          mean = paste0("Mean prob:", signif(mean, 2)*100, "%")
        ) %>%
        filter(Dataset == "UNITED" & Model == "Type 1") %>%
        mutate(
          Calculator = factor(Calculator, levels = c("Biomarkers", "No Biomarkers"), labels = c("Clinical features and biomarkers", "Clinical features"))
        ),
      mapping = aes(x = -Inf, y = -Inf, label = auc),
      size = 7,
      label.size = NA,
      hjust = -0.8,
      vjust = -0.5
    ) +
    theme(
      panel.spacing.x = unit(1.5, "lines")
    )),
  
  free(roc_curves %>%
    filter(Dataset == "UNITED" & Model == "Type 2") %>%
    ggplot(aes(x = 1- specificities, y = sensitivities)) +
    geom_path(
      data = roc_T2D_new_united_all_genes %>%
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
      mapping = aes(x = -Inf, y = -Inf, label = auc),
      size = 7,
      label.size = NA,
      hjust = -0.8,
      vjust = -0.5
    )),
  
  ncol = 1
  
) + patchwork::plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 11)
  )

pdf("figures/united_roc_thin_100.pdf", width = 8, height = 8)
plot_prob_rocs_united +
  patchwork::plot_layout(
    design = "
    AAAA
    #BB#
    "
  )
dev.off()




plot_prob_boxplot_united <- patchwork::wrap_plots(
  
  # Boxplots
  free(dataset.UNITED_type1_all_genes %>%
    select(M) %>%
    rename("Mody" = "M") %>%
    cbind(
      prob_with = colMeans(predictions_dataset.UNITED_type1_all_genes_with_T),
      prob_without = colMeans(predictions_dataset.UNITED_type1_all_genes_no_T)
    ) %>%
    gather("key", "Probability", -Mody) %>% 
    mutate(
      Mody = factor(Mody, levels = c(0, 1), labels = c("Negative", "Positive")),
      key = factor(key, levels = c("prob_without", "prob_with"), labels = c("Clinical features", "Clinical features and biomarkers"))
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
    )),
  
  # Boxplots
  free(dataset.UNITED_type2_all_genes %>%
    select(M) %>%
    rename("Mody" = "M") %>%
    cbind(
      Probability = colMeans(predictions_dataset.UNITED_type2_all_genes_new),
      key = ""
    ) %>%
    mutate(
      Mody = factor(Mody, levels = c(0, 1), labels = c("Negative", "Positive"))
    ) %>%
    ggplot() +
    geom_boxplot(aes(y = Probability, x = Mody)) +
    scale_y_continuous(labels = scales::percent) +
    theme_bw())
  
  , ncol = 1
  
) + patchwork::plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 11)
  )


pdf("figures/united_boxplot_thin_100.pdf", width = 8, height = 9)
plot_prob_boxplot_united +
  patchwork::plot_layout(
    design = "
    AAAA
    #BB#
    "
  )
dev.off()


################


prec_recal_curves <- data.frame(prob = colMeans(predictions_dataset.UNITED_type1_all_genes_with_T)) %>%
  cbind(Mody = dataset.UNITED_type1_all_genes$M) %>%
  pROC::roc(response = Mody, predictor = prob) %>%
  pROC::coords(ret = c("precision", "recall")) %>%
  as.data.frame() %>%
  mutate(
    ROCAUC = calc_auc_pr(dataset.UNITED_type1_all_genes$M, t(as.data.frame(colMeans(predictions_dataset.UNITED_type1_all_genes_with_T))), thinning = 1),
    mean = mean(colMeans(predictions_dataset.UNITED_type1_all_genes_with_T), na.rm = TRUE)
  ) %>%
  mutate(Dataset = "UNITED", Model = "Type 1", Calculator = "Biomarkers") %>%
  rbind(
    data.frame(prob = colMeans(predictions_dataset.UNITED_type1_all_genes_no_T)) %>%
      cbind(Mody = dataset.UNITED_type1_all_genes$M) %>%
      pROC::roc(response = Mody, predictor = prob) %>%
      pROC::coords(ret = c("precision", "recall")) %>%
      as.data.frame() %>%
      mutate(
        ROCAUC = calc_auc_pr(dataset.UNITED_type1_all_genes$M, t(as.data.frame(colMeans(predictions_dataset.UNITED_type1_all_genes_no_T))), thinning = 1),
        mean = mean(colMeans(predictions_dataset.UNITED_type1_all_genes_no_T), na.rm = TRUE)
      ) %>%
      mutate(Dataset = "UNITED", Model = "Type 1", Calculator = "No Biomarkers"), 
    data.frame(prob = colMeans(predictions_dataset.UNITED_type2_all_genes_new)) %>%
      cbind(Mody = dataset.UNITED_type2_all_genes$M) %>%
      pROC::roc(response = Mody, predictor = prob) %>%
      pROC::coords(ret = c("precision", "recall")) %>%
      as.data.frame() %>%
      mutate(
        ROCAUC = calc_auc_pr(dataset.UNITED_type2_all_genes$M, t(as.data.frame(colMeans(predictions_dataset.UNITED_type2_all_genes_new))), thinning = 1),
        mean = mean(colMeans(predictions_dataset.UNITED_type2_all_genes_new), na.rm = TRUE)
      ) %>%
      mutate(Dataset = "UNITED", Model = "Type 2", Calculator = " ")
  ) %>%
  rename("auc" = "ROCAUC")


dat_text <- prec_recal_curves %>%
  select(-precision, -recall) %>%
  distinct() %>%
  mutate(
    auc = paste0(" AUC:", signif(auc, 2), " "),
    mean = paste0("Mean prob:", signif(mean, 2)*100, "%"),
    Calculator = factor(Calculator, levels = c("Biomarkers", "No Biomarkers"), labels = c("Clinical features and biomarkers", "Clinical features"))
  )




plot_prob_prec_recal_united <- patchwork::wrap_plots(
  
  # roc insulin-treated
  free(prec_recal_curves %>%
         filter(Dataset == "UNITED" & Model == "Type 1") %>%
         mutate(
           Calculator = factor(Calculator, levels = c("No Biomarkers", "Biomarkers"), labels = c("Clinical features", "Clinical features and biomarkers")),
           iteration = 0
         ) %>%
         ggplot(aes(x = recall, y = precision)) +
         geom_path(
           data = prec_recal_T1D_no_T_united_all_genes %>%
             mutate(
               Calculator = "Clinical features"
             ) %>%
             rbind(
               prec_recal_T1D_with_T_united_all_genes %>%
                 mutate(Calculator = "Clinical features and biomarkers")
             ),
           aes(group = iteration), colour = "grey"
         ) + 
         geom_path() +
         theme_bw() +
         facet_grid(~factor(Calculator, levels = c("Clinical features", "Clinical features and biomarkers")), scales = "free",) +
         scale_y_continuous("Precision", labels = scales::percent) +
         scale_x_continuous("Recall", labels = scales::percent) +
         theme_bw() +
         geom_label(
           data = prec_recal_curves %>%
             select(-precision, -recall) %>%
             distinct() %>%
             mutate(
               auc = paste0(" AUC:", signif(auc, 2), " "),
               mean = paste0("Mean prob:", signif(mean, 2)*100, "%")
             ) %>%
             filter(Dataset == "UNITED" & Model == "Type 1") %>%
             mutate(
               Calculator = factor(Calculator, levels = c("Biomarkers", "No Biomarkers"), labels = c("Clinical features and biomarkers", "Clinical features"))
             ),
           mapping = aes(x = -Inf, y = -Inf, label = auc),
           size = 7,
           label.size = NA,
           hjust = -1,
           vjust = -6.5
         ) +
         theme(
           panel.spacing.x = unit(1.5, "lines")
         )),
  
  free(prec_recal_curves %>%
         filter(Dataset == "UNITED" & Model == "Type 2") %>%
         ggplot(aes(x = recall, y = precision)) +
         geom_path(
           data = prec_recal_T2D_new_united_all_genes %>%
             mutate(
               Calculator = "No Biomarkers"
             ),
           aes(group = iteration), colour = "grey"
         ) +
         geom_path() +
         theme_bw() +
         scale_y_continuous("Precision", labels = scales::percent) +
         scale_x_continuous("Recall", labels = scales::percent) +
         theme_bw() +
         geom_label(
           data = dat_text %>%
             filter(Dataset == "UNITED" & Model == "Type 2"),
           mapping = aes(x = -Inf, y = -Inf, label = auc),
           size = 7,
           label.size = NA,
           hjust = -0.9,
           vjust = -6.5
         )),
  
  ncol = 1
  
) + patchwork::plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 11)
  )

pdf("figures/united_prec_recal_thin_100.pdf", width = 8, height = 8)
plot_prob_prec_recal_united +
  patchwork::plot_layout(
    design = "
    AAAA
    #BB#
    "
  )
dev.off()



