#:--------------------------------------------------------
#   
# In this file we make predictions for UNITED insulin treated
#   using different antibody tests
#
#:--------------------------------------------------------

# load libraries
library(nimble)
library(rms)
library(tidyverse)
library(patchwork)

# load functions needed for generating data
source("data/create_data.R")

calculate_thresholds_diagnostics <- function(response, prediction) {
  
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
  
  ## select unique combinations of sensitivity and NTT (only the first occurance)
  matrix_thresholds <- matrix_thresholds %>%
    slice(which(duplicated(matrix_thresholds %>% select(-Thresholds, -`Pick-up rate`)) == FALSE))
  
  return(matrix_thresholds)
  
}



# Predictions
predictions_dataset.UNITED_type1_with_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_with_T.rds")
predictions_dataset.UNITED_type1_gad_with_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_gad_with_T.rds")
predictions_dataset.UNITED_type1_gad_ia2_with_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_gad_ia2_with_T.rds")




## Load population representative dataset


dataset.UNITED_type1 <- create_data(dataset = "united t1d", biomarkers = "reduced")

dataset.UNITED_type1_gad <- create_data(dataset = "united t1d", biomarkers = "full") %>%
  
  # check if the antibody variable in question is recorded
  mutate(A = GAD)


dataset.UNITED_type1_gad_ia2 <- create_data(dataset = "united t1d", biomarkers = "full") %>%
  
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


# Table information
thresholds_UNITED_t1d <- calculate_thresholds_diagnostics(dataset.UNITED_type1$M, predictions_dataset.UNITED_type1_with_T$prob)

thresholds_UNITED_t1d_gad <- calculate_thresholds_diagnostics(dataset.UNITED_type1_gad$M, predictions_dataset.UNITED_type1_gad_with_T$prob)

thresholds_UNITED_t1d_gad_ia2 <- calculate_thresholds_diagnostics(dataset.UNITED_type1_gad_ia2$M, predictions_dataset.UNITED_type1_gad_ia2_with_T$prob)



### 5%
thresholds_UNITED_t1d_gad %>%
  filter(Thresholds > 0.05) %>% head()

thresholds_UNITED_t1d_gad_ia2 %>%
  filter(Thresholds > 0.05) %>% head()

thresholds_UNITED_t1d %>%
  filter(Thresholds > 0.05) %>% head()

### 10%
thresholds_UNITED_t1d_gad %>%
  filter(Thresholds > 0.1) %>% head()

thresholds_UNITED_t1d_gad_ia2 %>%
  filter(Thresholds > 0.1) %>% head()

thresholds_UNITED_t1d %>%
  filter(Thresholds > 0.1) %>% head()

### 20%
thresholds_UNITED_t1d_gad %>%
  filter(Thresholds > 0.2) %>% head()

thresholds_UNITED_t1d_gad_ia2 %>%
  filter(Thresholds > 0.2) %>% head()

thresholds_UNITED_t1d %>%
  filter(Thresholds > 0.2) %>% head()

### 30%
thresholds_UNITED_t1d_gad %>%
  filter(Thresholds > 0.3) %>% head()

thresholds_UNITED_t1d_gad_ia2 %>%
  filter(Thresholds > 0.3) %>% head()

thresholds_UNITED_t1d %>%
  filter(Thresholds > 0.3) %>% head()




#################################

plot_antibody_boxplot <- patchwork::wrap_plots(
  
  # GAD
  dataset.UNITED_type1_gad %>%
    mutate(
      M = ifelse(is.na(M), 0, M)
    ) %>%
    select(A, M) %>%
    cbind(prob = predictions_dataset.UNITED_type1_gad_with_T$prob) %>%
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
    geom_boxplot(colour = c("white", "black", "black")) +
    geom_point(aes(alpha = M)) +
    scale_y_continuous("MODY probability", labels = scales::percent, breaks = seq(0, 0.7, 0.1), limits = c(0, 0.7)) +
    scale_alpha_manual(values = c(1, 0)) +
    ggtitle(paste0("GAD only tested (n=", nrow(dataset.UNITED_type1_gad %>% filter(!is.na(A))), ")")) +
    facet_wrap(~A_label) +
    theme_bw() +
    theme(
      legend.position = "none"
    ),
  
  # GAD & IA2
  dataset.UNITED_type1_gad_ia2 %>%
    mutate(
      M = ifelse(is.na(M), 0, M)
    ) %>%
    select(A, M) %>%
    cbind(prob = predictions_dataset.UNITED_type1_gad_ia2_with_T$prob) %>%
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
    geom_boxplot(colour = c("white", "black", "black")) +
    geom_point(aes(alpha = M)) +
    scale_y_continuous("MODY probability", labels = scales::percent, breaks = seq(0, 0.7, 0.1), limits = c(0, 0.7)) +
    scale_alpha_manual(values = c(1, 0)) +
    ggtitle(paste0("GAD & IA2 tested (n=", nrow(dataset.UNITED_type1_gad_ia2 %>% filter(!is.na(A))), ")")) +
    facet_wrap(~A_label) +
    theme_bw() +
    theme(
      legend.position = "none"
    ),
  
  # all antibodies
  dataset.UNITED_type1 %>%
    mutate(
      M = ifelse(is.na(M), 0, M)
    ) %>%
    select(A, M) %>%
    cbind(prob = predictions_dataset.UNITED_type1_with_T$prob) %>%
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
    geom_boxplot(colour = c("white", "black", "black")) +
    geom_point(aes(alpha = M)) +
    scale_y_continuous("MODY probability", labels = scales::percent, breaks = seq(0, 0.7, 0.1), limits = c(0, 0.7)) +
    scale_alpha_manual(values = c(1, 0)) +
    ggtitle(paste0("GAD & IA2 & ZnT8 tested (n=", nrow(dataset.UNITED_type1 %>% filter(!is.na(A))), ")")) +
    facet_wrap(~A_label) +
    theme_bw() +
    theme(
      legend.position = "none"
    ),
  
  nrow = 1
  
) +
  patchwork::plot_layout(axis_titles = "collect") &
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5)
  )


pdf("figures/plot_antibody_sensitivity.pdf", width = 13, height = 4)
plot_antibody_boxplot
dev.off()










