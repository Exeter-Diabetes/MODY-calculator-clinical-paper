##########################################################################################
#Figure 2


###############################################################################################
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
source("Data/create_data.R")
source("New_Data_Predictions/prediction_functions.R")

# load files required
predictions_UNITED_type1_no_T_full <- 
  readRDS("Model_Predictions/predictions_dataset.UNITED_type1_all_genes_no_T_full.rds")
predictions_UNITED_type1_with_T_full  <- 
  readRDS("Model_Predictions/predictions_dataset.UNITED_type1_all_genes_with_T_full.rds")
predictions_UNITED_type2_full <- 
  readRDS("Model_Predictions/predictions_dataset.UNITED_type2_all_genes_new_full.rds")

predictions_UNITED_type1_no_T <- 
  readRDS("Model_Predictions/predictions_dataset.UNITED_type1_all_genes_no_T.rds")
predictions_UNITED_type1_with_T <- 
  readRDS("Model_Predictions/predictions_dataset.UNITED_type1_all_genes_with_T.rds")
predictions_UNITED_type2 <- 
  readRDS("Model_Predictions/predictions_dataset.UNITED_type2_all_genes_new.rds")

# load datasets
## Load population representative dataset
dataset_UNITED_type1 <- create_data(dataset = "united t1d", 
                                              commonmody = FALSE) %>%
  
  ## if MODY testing missing, change to 0
  mutate(M = ifelse(is.na(M), 0, M))

dataset_UNITED_type2 <- create_data(dataset = "united t2d", 
                                              commonmody = FALSE)

#merge probs to datasets
dataset_UNITED_type1 <- cbind(dataset_UNITED_type1, predictions_UNITED_type1_with_T)
dataset_UNITED_type2 <- cbind(dataset_UNITED_type2, predictions_UNITED_type2)

#merge to joint dataset
UNITED_joint <- full_join(dataset_UNITED_type1, dataset_UNITED_type2)



#Figure prep --------------------------------------------------------------------------------------

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
auc_UNITED_type1_no_T <- calc_auroc(dataset_UNITED_type1$M, 
                                      predictions_UNITED_type1_no_T_full,
                                      thinning = 10)
# quantile(auc_UNITED_type1_no_T, probs = c(0.025, 0.5, 0.975)) # thinning = 10
# 2.5%       50%     97.5% 
# 0.7053243 0.7924399 0.8190743 

### Biomarker models
auc_UNITED_type1_with_T <- calc_auroc(dataset_UNITED_type1$M, 
                                      predictions_UNITED_type1_with_T_full,
                                      thinning = 10)
# quantile(auc_UNITED_type1_with_T, probs = c(0.025, 0.5, 0.975)) # thinning = 10
# 2.5%       50%     97.5% 
# 0.9487113 0.9768041 0.9779210 

## Type 2 UNITED
auc_UNITED_type2 <- calc_auroc(dataset_UNITED_type2$M,
                               predictions_UNITED_type2_full,
                               thinning = 10)
# quantile(auc_UNITED_type2, probs = c(0.025, 0.5, 0.975)) # thinning = 10
# 2.5%       50%     97.5% 
# 0.8465473 0.8618926 0.8768116  

#:------------------------------------------------------------
roc_curves <- data.frame(prob = colMeans(predictions_UNITED_type1_with_T_full)) %>%
  cbind(Mody = dataset_UNITED_type1$M) %>%
  pROC::roc(response = Mody, predictor = prob) %>%
  magrittr::extract(2:3) %>%
  as.data.frame() %>%
  mutate(
    auc =  unname(data.frame(prob = colMeans(predictions_UNITED_type1_with_T_full)) %>%
                    cbind(Mody = dataset_UNITED_type1$M) %>%
                    pROC::roc(response = Mody, predictor = prob) %>%
                    magrittr::extract(c(9)) %>%
                    unlist()),
    auc_low = quantile(auc_UNITED_type1_with_T, probs = c(0.025)),
    auc_high = quantile(auc_UNITED_type1_with_T, probs= c(0.975)),
    mean = mean(colMeans(predictions_UNITED_type1_with_T_full), na.rm = TRUE)
  ) %>%
  mutate(Dataset = "UNITED", Model = "Type 1", Calculator = "Biomarkers") %>%
  rbind(
    data.frame(prob = colMeans(predictions_UNITED_type1_no_T_full)) %>%
      cbind(Mody = dataset_UNITED_type1$M) %>%
      pROC::roc(response = Mody, predictor = prob) %>%
      magrittr::extract(2:3) %>%
      as.data.frame() %>%
      mutate(
        auc = unname(data.frame(prob = colMeans(predictions_UNITED_type1_no_T_full)) %>%
                       cbind(Mody = dataset_UNITED_type1$M) %>%
                       pROC::roc(response = Mody, predictor = prob) %>%
                       magrittr::extract(c(9)) %>%
                       unlist()),
        auc_low = quantile(auc_UNITED_type1_no_T, probs = c(0.025)),
        auc_high = quantile(auc_UNITED_type1_no_T, probs= c(0.975)),
        mean = mean(colMeans(predictions_UNITED_type1_no_T_full), na.rm = TRUE)
      ) %>%
      mutate(Dataset = "UNITED", Model = "Type 1", Calculator = "No Biomarkers"), 
    data.frame(prob = colMeans(predictions_UNITED_type2_full)) %>%
      cbind(Mody = dataset_UNITED_type2$M) %>%
      pROC::roc(response = Mody, predictor = prob) %>%
      magrittr::extract(2:3) %>%
      as.data.frame() %>%
      mutate(
        auc = unname(data.frame(prob = colMeans(predictions_UNITED_type2_full)) %>%
                       cbind(Mody = dataset_UNITED_type2$M) %>%
                       pROC::roc(response = Mody, predictor = prob) %>%
                       magrittr::extract(c(9)) %>%
                       unlist()),
        auc_low = quantile(auc_UNITED_type2, probs = c(0.025)),
        auc_high = quantile(auc_UNITED_type2, probs= c(0.975)),
        mean = mean(colMeans(predictions_UNITED_type2_full), na.rm = TRUE)
      ) %>%
      mutate(Dataset = "UNITED", Model = "Type 2", Calculator = "No Biomarkers")
  )


dat_text <- roc_curves %>%
  select(-sensitivities, -specificities) %>%
  distinct() %>%
  mutate(
    auc_full = paste0("AUC: ", signif(auc, 2), " [", signif(auc_low, 2), "-", signif(auc_high, 2), "]"),
    mean = paste0("Mean prob:", signif(mean, 2)*100, "%"),
    Calculator = factor(Calculator, 
                        levels = c("Biomarkers", 
                                   "No Biomarkers"), 
                        labels = c("Clinical features and biomarkers", 
                                   "Clinical features"))
  )



# Plots ------------------------------------------------------------------------------
plot_prob_fig2 <- patchwork::wrap_plots(
  # Row 1 - early insulin-treated
  patchwork::wrap_plots(
    #A.1 - ROC of clinicial feautres & no biomarkers
    roc_curves %>%
      filter(Dataset == "UNITED" & 
               Model == "Type 1" & 
               Calculator == "No Biomarkers") %>%
      mutate(
        Calculator = factor(Calculator, 
                            levels = c("No Biomarkers"), 
                            labels = c("Clinical features")),
        iteration = 0
      ) %>%
      ggplot(aes(x = 1- specificities, y = sensitivities)) +
      geom_path() +
      theme_bw() +
      facet_grid(~factor(Calculator, 
                         levels = c("Clinical features"), 
                         labels = c("Early-insulin-treated: clinical features, no biomarkers")), 
                 scales = "free",
                 labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
      scale_y_continuous("Sensitivity", labels = scales::percent) +
      scale_x_continuous("1- Specificity", labels = scales::percent) +
      theme_bw() +
      geom_label(
        data = dat_text %>%
          filter(Dataset == "UNITED" & 
                   Model == "Type 1" & 
                   Calculator == "Clinical features"),
        mapping = aes(x = 0.55, y = 0.1, label = auc_full, hjust = "center"),
        size = 7,
        label.r = unit(0, "pt"),
        label.padding=unit(0.4, "lines")
      ) +
      theme(
        panel.spacing.x = unit(1.5, "lines")
      ),
    #A.2 - ROC of clinical features & biomarkers
    roc_curves %>%
      filter(Dataset == "UNITED" & 
               Model == "Type 1" & 
               Calculator == "Biomarkers") %>%
      mutate(
        Calculator = factor(Calculator, 
                            levels = c("Biomarkers"), 
                            labels = c("Clinical features and biomarkers")),
        iteration = 0
      ) %>%
      ggplot(aes(x = 1- specificities, y = sensitivities)) +
      geom_path() +
      theme_bw() +
      facet_grid(~factor(Calculator, 
                         levels = c("Clinical features and biomarkers"), 
                         labels = c("Early-insulin-treated: clinical features, with biomarkers")), 
                 scales = "free",labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
      scale_y_continuous("Sensitivity", labels = scales::percent) +
      scale_x_continuous("1- Specificity", labels = scales::percent) +
      theme_bw() +
      geom_label(
        data = dat_text %>%
          filter(Dataset == "UNITED" & 
                   Model == "Type 1" & 
                   Calculator == "Clinical features and biomarkers"),
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
  #Row 2 - not early insulin treated ROC (B) & distribution (C)
  patchwork::wrap_plots(
    #B - ROC of not-early-insulin-treated
    roc_curves %>%
      filter(Dataset == "UNITED" & Model == "Type 2") %>%
      mutate(Calculator = factor(Calculator, 
                                 levels = "No Biomarkers", 
                                 labels = "Not-early-insulin-treated")) %>%
      ggplot(aes(x = 1- specificities, y = sensitivities)) +
      geom_path() +
      theme_bw() +
      facet_grid(~factor(Calculator, 
                         levels = c("Not-early-insulin-treated"), 
                         labels = c("Not-early-insulin-treated: clinical features, no biomarkers")), 
                 scales = "free",
                 labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
      scale_y_continuous("Sensitivity", labels = scales::percent) +
      scale_x_continuous("1- Specificity", labels = scales::percent) +
      theme_bw() +
      geom_label(
        data = dat_text %>%
          filter(Dataset == "UNITED" & Model == "Type 2") %>%
          mutate(
            Calculator = factor(Calculator, 
                                levels = c("Clinical features"), 
                                labels = c("Not-early-insulin-treated"))
          ),
        mapping = aes(x = 0.55, y = 0.1, label = auc_full, hjust = "center"),
        size = 7,
        label.r = unit(0, "pt"),
        label.padding=unit(0.4, "lines")
      ),
    #C - Distribution curves for joint dataset
    patchwork::wrap_plots(
      #Density plot
      UNITED_joint %>%
        mutate(prob = UNITED_joint$prob) %>%
        filter(M == 0) %>%
        select(M, prob) %>%
        rename("Mody" = "M") %>%
        mutate(
          Mody = factor(Mody, levels = c(0, 1), labels = c("Non-MODY", "MODY")),
        ) %>%
        ggplot() +
        geom_density(aes(x = prob), fill = "grey") +
        geom_vline(xintercept = 0.05) +
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
        ylab("Non-MODY (n=1266)"),
      #point
      UNITED_joint %>%
        mutate(prob = UNITED_joint$prob) %>%
        filter(M == 1) %>%
        select(M, prob) %>%
        rename("Mody" = "M") %>%
        mutate(
          Mody = factor(Mody, levels = c(0, 1), labels = c("Non-MODY", "MODY")),
        ) %>%
        ggplot() +
        geom_point(aes(x=prob, y=0), position = position_jitter(height = 0.1, seed = 20)) +
        geom_vline(xintercept = 0.05) +
        coord_cartesian(xlim =c(0, 1), ylim = c(-0.15, 0.15)) +
        scale_x_continuous(labels = scales::percent) +
        theme_classic() +
        theme(
          axis.ticks.y = element_blank(), 
          axis.text.y = element_blank(), 
          #axis.title.y = element_blank()
        ) +
        ylab("MODY \n cases \n (n=33)") +
        xlab("Model probabilities"),
      ncol = 1, nrow = 2, heights = c(4,1)
  ),
  ncol = 2, nrow = 1
  ),
  ncol = 1
  ) +

patchwork::plot_annotation(tag_levels = list(c("A.1", "A.2", "B", "C"))) &
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 11)
  )

# Making plots
pdf("Figures/Figure2.pdf", width = 13, height = 9)
plot_prob_fig2
dev.off()


# Making plots for presenting
ggsave("Figures/Figure2.tif", width = 13, height = 9, dpi= 1000)
plot_prob_fig2
dev.off()