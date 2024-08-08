#:--------------------------------------------------------
#   
# In this file we make predictions for UNITED insulin treated
#   using different antibody tests
#
#:--------------------------------------------------------

# load libraries
library(nimble)
library(rms)
library(pROC)
library(tidyverse)
library(patchwork)

# load functions needed for generating data
source("data/create_data.R")

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



# Predictions
predictions_dataset.UNITED_type1_with_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_with_T.rds")
predictions_dataset.UNITED_type1_with_T_full <- readRDS("model_predictions/predictions_dataset.UNITED_type1_with_T_full.rds")
predictions_dataset.UNITED_type1_all_genes_with_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_all_genes_with_T.rds")
predictions_dataset.UNITED_type1_all_genes_with_T_full <- readRDS("model_predictions/predictions_dataset.UNITED_type1_all_genes_with_T_full.rds")
predictions_dataset.UNITED_type1_gad_all_genes_with_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_gad_all_genes_with_T.rds")
predictions_dataset.UNITED_type1_gad_all_genes_with_T_full <- readRDS("model_predictions/predictions_dataset.UNITED_type1_gad_all_genes_with_T_full.rds")
predictions_dataset.UNITED_type1_gad_ia2_all_genes_with_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_gad_ia2_all_genes_with_T.rds")
predictions_dataset.UNITED_type1_gad_ia2_all_genes_with_T_full <- readRDS("model_predictions/predictions_dataset.UNITED_type1_gad_ia2_all_genes_with_T_full.rds")
predictions_dataset.UNITED_type1_no_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_no_T.rds")
predictions_dataset.UNITED_type1_sensitivity_analysis_with_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_sensitivity_analysis_with_T.rds")


## Load population representative dataset


dataset.UNITED_type1 <- create_data(dataset = "united t1d", biomarkers = "reduced")

dataset.UNITED_type1_all_genes <- create_data(dataset = "united t1d", biomarkers = "reduced", commonmody = FALSE)

dataset.UNITED_type1_gad_all_genes <- create_data(dataset = "united t1d", biomarkers = "full", commonmody = FALSE) %>%
  
  # check if the antibody variable in question is recorded
  mutate(A = GAD)


dataset.UNITED_type1_gad_ia2_all_genes <- create_data(dataset = "united t1d", biomarkers = "full", commonmody = FALSE) %>%
  
  # check if the antibody variable in question is recorded
  mutate(
    A = ifelse(is.na(GAD) & is.na(IA2), NA,
               ifelse(!is.na(GAD) & (GAD == 1 & is.na(IA2)), 1,
                      ifelse(!is.na(GAD) & (GAD == 0 & is.na(IA2)), 0,
                             ifelse(!is.na(IA2) & (IA2 == 1 & is.na(GAD)), 1,
                                    ifelse(!is.na(IA2) & (IA2 == 0 & is.na(GAD)), 0,
                                           ifelse(IA2 == 1 | GAD == 1, 1, 0))))))
  )





####################################################################################################################################################################################################################
####################################################################################################################################################################################################################

####################################################################################################################################################################################################################
####################################################################################################################################################################################################################

# Sensitivity analysis: model fitted with GAD information only vs with all antibody information

plot_differences <- patchwork::wrap_plots(
  data.frame(
  prob_sens = predictions_dataset.UNITED_type1_sensitivity_analysis_with_T$prob,
  prob = predictions_dataset.UNITED_type1_with_T$prob
) %>%
  ggplot(aes(x = prob, y = prob_sens)) +
  geom_point() +
  xlab("All information probabilities") +
  ylab("GAD only information probabilities") +
  geom_abline(aes(intercept = 0, slope = 1)) +
  theme_bw(),

  data.frame(
    diff = predictions_dataset.UNITED_type1_sensitivity_analysis_with_T$prob - predictions_dataset.UNITED_type1_with_T$prob
  ) %>%
  ggplot(aes(x = diff)) +
  geom_histogram() +
  xlab("GAD only probabilities - All information probabilities") +
  theme_bw()

)






####################################################################################################################################################################################################################
####################################################################################################################################################################################################################

####################################################################################################################################################################################################################
####################################################################################################################################################################################################################


# Table information
thresholds_UNITED_t1d_no_T <- calculate_thresholds_diagnostics(dataset.UNITED_type1$M, predictions_dataset.UNITED_type1_no_T$prob)

thresholds_UNITED_t1d_all_genes <- calculate_thresholds_diagnostics(dataset.UNITED_type1_all_genes$M, predictions_dataset.UNITED_type1_all_genes_with_T$prob)

thresholds_UNITED_t1d_gad_all_genes <- calculate_thresholds_diagnostics(dataset.UNITED_type1_gad_all_genes$M, predictions_dataset.UNITED_type1_gad_all_genes_with_T$prob)

thresholds_UNITED_t1d_gad_ia2_all_genes <- calculate_thresholds_diagnostics(dataset.UNITED_type1_gad_ia2_all_genes$M, predictions_dataset.UNITED_type1_gad_ia2_all_genes_with_T$prob)


### 5%
thresholds_UNITED_t1d_no_T %>%
  filter(Thresholds == 0.05) %>% arrange(Thresholds) %>% head()

thresholds_UNITED_t1d_gad_all_genes %>%
  filter(Thresholds == 0.05) %>% arrange(Thresholds) %>% head()

thresholds_UNITED_t1d_gad_ia2_all_genes %>%
  filter(Thresholds == 0.05) %>% arrange(Thresholds)  %>% head()

thresholds_UNITED_t1d_all_genes %>%
  filter(Thresholds == 0.05) %>% arrange(Thresholds)  %>% head()

### 10%
thresholds_UNITED_t1d_no_T %>%
  filter(Thresholds == 0.1) %>% arrange(Thresholds) %>% head()

thresholds_UNITED_t1d_gad_all_genes %>%
  filter(Thresholds == 0.1) %>% arrange(Thresholds)  %>% head()

thresholds_UNITED_t1d_gad_ia2_all_genes %>%
  filter(Thresholds == 0.1) %>% arrange(Thresholds)  %>% head()

thresholds_UNITED_t1d_all_genes %>%
  filter(Thresholds == 0.1) %>% arrange(Thresholds)  %>% head()

### 20%
thresholds_UNITED_t1d_no_T %>%
  filter(Thresholds == 0.2) %>% arrange(Thresholds) %>% head()

thresholds_UNITED_t1d_gad_all_genes %>%
  filter(Thresholds == 0.2) %>% arrange(Thresholds)  %>% head()

thresholds_UNITED_t1d_gad_ia2_all_genes %>%
  filter(Thresholds == 0.2) %>% arrange(Thresholds)  %>% head()

thresholds_UNITED_t1d_all_genes %>%
  filter(Thresholds == 0.2) %>% arrange(Thresholds)  %>% head()

### 30%
thresholds_UNITED_t1d_no_T %>%
  filter(Thresholds == 0.3) %>% arrange(Thresholds) %>% head()

thresholds_UNITED_t1d_gad_all_genes %>%
  filter(Thresholds == 0.3) %>% arrange(Thresholds)  %>% head()

thresholds_UNITED_t1d_gad_ia2_all_genes %>%
  filter(Thresholds == 0.3) %>% arrange(Thresholds)  %>% head()

thresholds_UNITED_t1d_all_genes %>%
  filter(Thresholds == 0.3) %>% arrange(Thresholds)  %>% head()


#################################

# Roc curves

### GAD only 
roc_T1D_with_T_united_gad_all_genes <- calc_roc(dataset.UNITED_type1_gad_all_genes %>% mutate(M = ifelse(is.na(M), 0, M)) %>% select(M) %>% unlist(), predictions_dataset.UNITED_type1_gad_all_genes_with_T_full, thinning = 1000)

# AUC ROC
auc_roc_T1D_with_T_united_gad_all_genes <- unname(data.frame(prob = predictions_dataset.UNITED_type1_gad_all_genes_with_T$prob) %>%
                                                    cbind(Mody = dataset.UNITED_type1_gad_all_genes$M) %>%
                                                    mutate(Mody = ifelse(is.na(Mody), 0, Mody)) %>%
                                                    pROC::roc(response = Mody, predictor = prob) %>%
                                                    magrittr::extract(c(9)) %>%
                                                    unlist())

# plot for ROC with grey being iterations, black being the ROC for average prediction
plot_roc_T1D_with_T_united_gad_all_genes <- ggplot() +
  ## all iterations
  geom_path(
    data = roc_T1D_with_T_united_gad_all_genes,
    aes(x = 1-specificities, y= sensitivities, group = iteration), colour = "grey"
  ) +
  ## average predictions
  geom_path(
    data = pROC::roc(response = dataset.UNITED_type1_gad_all_genes %>% mutate(M = ifelse(is.na(M), 0, M)) %>% select(M) %>% unlist(), predictor = predictions_dataset.UNITED_type1_gad_all_genes_with_T$prob, levels = c(0,1), direction = "<") %>%
      magrittr::extract(c(2:3)) %>%
      as.data.frame(),
    aes(x = 1-specificities, y= sensitivities), colour = "black"
  ) +
  theme_bw() +
  scale_y_continuous("Sensitivity", labels = scales::percent) +
  scale_x_continuous("1- Specificity", labels = scales::percent)


### GAD and IA2
roc_T1D_with_T_united_gad_ia2_all_genes <- calc_roc(dataset.UNITED_type1_gad_ia2_all_genes %>% mutate(M = ifelse(is.na(M), 0, M)) %>% select(M) %>% unlist(), predictions_dataset.UNITED_type1_gad_ia2_all_genes_with_T_full, thinning = 1000)

# AUC ROC
auc_roc_T1D_with_T_united_gad_ia2_all_genes <- unname(data.frame(prob = predictions_dataset.UNITED_type1_gad_ia2_all_genes_with_T$prob) %>%
                                                    cbind(Mody = dataset.UNITED_type1_gad_ia2_all_genes$M) %>%
                                                    mutate(Mody = ifelse(is.na(Mody), 0, Mody)) %>%
                                                    pROC::roc(response = Mody, predictor = prob) %>%
                                                    magrittr::extract(c(9)) %>%
                                                    unlist())

# plot for ROC with grey being iterations, black being the ROC for average prediction
plot_roc_T1D_with_T_united_gad_ia2_all_genes <- ggplot() +
  ## all iterations
  geom_path(
    data = roc_T1D_with_T_united_gad_ia2_all_genes,
    aes(x = 1-specificities, y= sensitivities, group = iteration), colour = "grey"
  ) +
  ## average predictions
  geom_path(
    data = pROC::roc(response = dataset.UNITED_type1_gad_ia2_all_genes %>% mutate(M = ifelse(is.na(M), 0, M)) %>% select(M) %>% unlist(), predictor = predictions_dataset.UNITED_type1_gad_ia2_all_genes_with_T$prob, levels = c(0,1), direction = "<") %>%
      magrittr::extract(c(2:3)) %>%
      as.data.frame(),
    aes(x = 1-specificities, y= sensitivities), colour = "black"
  ) +
  theme_bw() +
  scale_y_continuous("Sensitivity", labels = scales::percent) +
  scale_x_continuous("1- Specificity", labels = scales::percent)

### All antibodies
roc_T1D_with_T_united_all_genes <- calc_roc(dataset.UNITED_type1_all_genes %>% mutate(M = ifelse(is.na(M), 0, M)) %>% select(M) %>% unlist(), predictions_dataset.UNITED_type1_all_genes_with_T_full, thinning = 1000)

# AUC ROC
auc_roc_T1D_with_T_united_all_genes <- unname(data.frame(prob = predictions_dataset.UNITED_type1_all_genes_with_T$prob) %>%
                                                        cbind(Mody = dataset.UNITED_type1_all_genes$M) %>%
                                                        mutate(Mody = ifelse(is.na(Mody), 0, Mody)) %>%
                                                        pROC::roc(response = Mody, predictor = prob) %>%
                                                        magrittr::extract(c(9)) %>%
                                                        unlist())

# plot for ROC with grey being iterations, black being the ROC for average prediction
plot_roc_T1D_with_T_united_all_genes <- ggplot() +
  ## all iterations
  geom_path(
    data = roc_T1D_with_T_united_all_genes,
    aes(x = 1-specificities, y= sensitivities, group = iteration), colour = "grey"
  ) +
  ## average predictions
  geom_path(
    data = pROC::roc(response = dataset.UNITED_type1_all_genes %>% mutate(M = ifelse(is.na(M), 0, M)) %>% select(M) %>% unlist(), predictor = predictions_dataset.UNITED_type1_all_genes_with_T$prob, levels = c(0,1), direction = "<") %>%
      magrittr::extract(c(2:3)) %>%
      as.data.frame(),
    aes(x = 1-specificities, y= sensitivities), colour = "black"
  ) +
  theme_bw() +
  scale_y_continuous("Sensitivity", labels = scales::percent) +
  scale_x_continuous("1- Specificity", labels = scales::percent)


#################################

plot_antibody_boxplot <- patchwork::wrap_plots(
  
  # GAD
  ## boxplot
  dataset.UNITED_type1_gad_all_genes %>%
    mutate(
      M = ifelse(is.na(M), 0, M)
    ) %>%
    select(A, M) %>%
    cbind(prob = predictions_dataset.UNITED_type1_gad_all_genes_with_T$prob) %>%
    drop_na() %>%
    mutate(
      total = n(),
      M = ifelse(M == 0, "Non-MODY", "MODY")
    ) %>%
    group_by(A) %>%
    mutate(
      A_count = n(),
      A = ifelse(A == 0, "Negative", "Positive"),
      A_label = paste0(A, "\nn=", A_count, " (", plyr::round_any(A_count/total*100, accuracy = 1, f = round),"%)"),
      point_alpha = ifelse(M == 1, 1, 0)
    ) %>%
    ungroup() %>%
    select(-A, -total, -A_count) %>%
    ggplot(aes(x = M, y = prob)) +
    geom_boxplot(colour = c("white", "black", "black"), alpha = c(0, 1, 1)) +
    geom_point(aes(alpha = M)) +
    scale_y_continuous("MODY probability", labels = scales::percent, breaks = seq(0, 0.7, 0.1), limits = c(0, 0.7)) +
    scale_alpha_manual(values = c(1, 0)) +
    ggtitle(paste0("GAD only tested (n=", nrow(dataset.UNITED_type1_gad_all_genes %>% filter(!is.na(A))), ")")) +
    facet_wrap(~A_label) +
    theme_bw() +
    theme(
      legend.position = "none"
    ),
  
  ## roc curve
  plot_roc_T1D_with_T_united_gad_all_genes +
    geom_label(
      mapping = aes(x = -Inf, y = -Inf), label = paste0("AUC:", signif(auc_roc_T1D_with_T_united_gad_all_genes, 2)),
      size = 7,
      label.size = NA,
      hjust = -1,
      vjust = -0.5
    ),
  
  # GAD & IA2
  dataset.UNITED_type1_gad_ia2_all_genes %>%
    mutate(
      M = ifelse(is.na(M), 0, M)
    ) %>%
    select(A, M) %>%
    cbind(prob = predictions_dataset.UNITED_type1_gad_ia2_all_genes_with_T$prob) %>%
    drop_na() %>%
    mutate(
      total = n(),
      M = ifelse(M == 0, "Non-MODY", "MODY")
    ) %>%
    group_by(A) %>%
    mutate(
      A_count = n(),
      A = ifelse(A == 0, "Negative", "Positive"),
      A_label = paste0(A, "\nn=", A_count, " (", plyr::round_any(A_count/total*100, accuracy = 1, f = round),"%)"),
      point_alpha = ifelse(M == 1, 1, 0)
    ) %>%
    ungroup() %>%
    select(-A, -total, -A_count) %>%
    ggplot(aes(x = M, y = prob)) +
    geom_boxplot(colour = c("white", "black", "black"), alpha = c(0, 1, 1)) +
    geom_point(aes(alpha = M)) +
    scale_y_continuous("MODY probability", labels = scales::percent, breaks = seq(0, 0.7, 0.1), limits = c(0, 0.7)) +
    scale_alpha_manual(values = c(1, 0)) +
    ggtitle(paste0("GAD & IA2 tested (n=", nrow(dataset.UNITED_type1_gad_ia2_all_genes %>% filter(!is.na(A))), ")")) +
    facet_wrap(~A_label) +
    theme_bw() +
    theme(
      legend.position = "none"
    ),
  
  ## roc curve
  plot_roc_T1D_with_T_united_gad_ia2_all_genes +
    geom_label(
      mapping = aes(x = -Inf, y = -Inf), label = paste0("AUC:", signif(auc_roc_T1D_with_T_united_gad_ia2_all_genes, 2)),
      size = 7,
      label.size = NA,
      hjust = -1,
      vjust = -0.5
    ),
  
  # all antibodies
  dataset.UNITED_type1_all_genes %>%
    mutate(
      M = ifelse(is.na(M), 0, M)
    ) %>%
    select(A, M) %>%
    cbind(prob = predictions_dataset.UNITED_type1_all_genes_with_T$prob) %>%
    drop_na() %>%
    mutate(
      total = n(),
      M = ifelse(M == 0, "Non-MODY", "MODY")
    ) %>%
    group_by(A) %>%
    mutate(
      A_count = n(),
      A = ifelse(A == 0, "Negative", "Positive"),
      A_label = paste0(A, "\nn=", A_count, " (", plyr::round_any(A_count/total*100, accuracy = 1, f = round),"%)"),
      point_alpha = ifelse(M == 1, 1, 0)
    ) %>%
    ungroup() %>%
    select(-A, -total, -A_count) %>%
    ggplot(aes(x = M, y = prob)) +
    geom_boxplot(colour = c("white", "black", "black"), alpha = c(0, 1, 1)) +
    geom_point(aes(alpha = M)) +
    scale_y_continuous("MODY probability", labels = scales::percent, breaks = seq(0, 0.7, 0.1), limits = c(0, 0.7)) +
    scale_alpha_manual(values = c(1, 0)) +
    ggtitle(paste0("GAD & IA2 & ZnT8 tested (n=", nrow(dataset.UNITED_type1_all_genes %>% filter(!is.na(A))), ")")) +
    facet_wrap(~A_label) +
    theme_bw() +
    theme(
      legend.position = "none"
    ),
  
  ## roc curve
  plot_roc_T1D_with_T_united_all_genes +
    geom_label(
      mapping = aes(x = -Inf, y = -Inf), label = paste0("AUC:", signif(auc_roc_T1D_with_T_united_all_genes, 2)),
      size = 7,
      label.size = NA,
      hjust = -1,
      vjust = -0.5
    ),
  
  nrow = 3
  
) +
  patchwork::plot_annotation(tag_levels = list(c("A.1", "A.2", "B.1", "B.2", "C.1", "C.2"))) +
  patchwork::plot_layout(axis_titles = "collect") &
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5)
  )


pdf("figures/plot_antibody_sensitivity.pdf", width = 9, height = 12)
plot_antibody_boxplot
dev.off()










