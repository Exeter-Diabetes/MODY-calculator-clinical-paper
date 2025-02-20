##########################################################################################
#Figure 3


###############################################################################################

# load libraries
library(tidyverse)
library(nimble)
library(pROC)
library(PRROC)

# load functions needed
source("data/create_data.R")
source("new_data_predictions/prediction_functions.R")

#readin datasets --------------------------------------------------------------------------------------
##data preparation ----------------------------------------------------------------------------------------------------------------------
### read in data ------------------------------------------------------------------------------------------------------------------------
#readin MYDIABETES
#load T1D
load("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper/MY_T1D.RData")
#load T2D 
load("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper/MY_T2D.RData")

# MyDiabetes ------------------------------------------------------------------------------------
#load T1D
load("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper/MY_T1D.RData")
#load T2D 
load("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper/MY_T2D.RData")
# LIMIT TO ONLY WHITE ETHNICITY
MY_T1D <- MY_T1D %>%
  filter(ethnicity == "White")
MY_T1D <- MY_T1D %>%
  rename(
    hba1c_mmol = hba1c
  )

MY_T1D <- MY_T1D %>%
  rename(
    hba1c = hba1c_perc, 
    bmi = BMI
  )

#T2D data
MY_T2D <- MY_T2D %>%
  filter(ethnicity == "White")

MY_T2D <- MY_T2D %>%
  rename(
    hba1c_mmol = hba1c
  )

MY_T2D <- MY_T2D %>%
  rename(
    hba1c = hba1c_perc, 
    bmi = BMI
  )

# UNITED paediatric ---------------------------------------------------------------------
dataset.UNITED_p <- create_data(dataset = "united t1d pediatrics", commonmody = FALSE)

#need to change M=NA to M=0
dataset.UNITED_p <- dataset.UNITED_p %>%
  mutate(M = ifelse(is.na(M) == TRUE, 0, M)) %>%
  filter(!is.na(T))
#checked if worked: should have M=1 (n=7) & M=0 (n=1164)
table(dataset.UNITED_p$M)

#T2D
dataset.UNITED_type2p <- dataset.UNITED_p %>%
  filter(tti == "" | tti == "Greater than 12 months")
table(dataset.UNITED_type2p$M)

#T1D
dataset.UNITED_type1p <- dataset.UNITED_p %>%
  filter(!(id %in% dataset.UNITED_type2p$id))
table(dataset.UNITED_type1p$M)

# External dataset -----------------------------------------------------------------------

# Early-insulin-treated
dataset.MYDIABETES_type1 <- MY_T1D %>%
  mutate(study = "MYDIABETES") %>%
  select(MY_ID, study, agerec, agedx, sex, bmi, pardm, insoroha, hba1c, C, A, M, Gene, biomark_status) 

dataset.UNITED_type1p <- dataset.UNITED_type1p %>%
  mutate(study = "UNITED paediatric") %>%
  filter(!is.na(T))
dataset_type1 <- full_join(dataset.MYDIABETES_type1, dataset.UNITED_type1p, by = c("study","agerec", "agedx", "sex", "bmi", "pardm", "insoroha", "hba1c", "C", "A", "M", "Gene")) %>%
  mutate(M = ifelse(is.na(M), 0, M), 
         Gene = ifelse(Gene == "", NA, Gene), 
         biomark_status = ifelse(is.na(biomark_status), 
                                 ifelse(C == 1 & A == 0, 0, 1),
                                 biomark_status))
dataset_type1 <- dataset_type1 %>%
  rename(agerec_original = agerec) %>%
  mutate(agerec = ifelse(agerec_original > 50, 50, agerec_original))
# Not-early-insulin-treated
dataset.MYDIABETES_type2 <- MY_T2D %>%
  mutate(study = "MYDIABETES") %>%
  select(MY_ID, study, agerec, agedx, sex, bmi, pardm, insoroha, hba1c, C, A, M, Gene, biomark_status) 

dataset.UNITED_type2p <- dataset.UNITED_type2p %>%
  mutate(study = "UNITED paediatric") %>%
  select(id, study, agerec, agedx, sex, bmi, pardm, insoroha, hba1c, C, A, M, Gene) 

dataset_type2 <- full_join(dataset.MYDIABETES_type2, dataset.UNITED_type2p, by = c("study","agerec", "agedx", "sex", "bmi", "pardm", "insoroha", "hba1c", "C", "A", "M", "Gene")) %>%
  mutate(M = ifelse(is.na(M), 0, M), 
         Gene = ifelse(Gene == "", NA, Gene), 
         biomark_status = ifelse(is.na(biomark_status), 
                                 ifelse(C == 1 & A == 0, 0, 1),
                                 biomark_status)) 



#run predictions ------------------------------------------------------------------------------
## In MYDIABETES
setwd("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper/new_data_predictions")
# ## load posteriors
rcs_parms <- readRDS("rcs_parms.rds")

##T1D/EARLY-INSULIN-TREATED MODEL -------------------------------------------------------------------
## IN UNITED PAEDIATRICS
posterior_samples_T1D <- readRDS("type_1_model_posteriors_thin_100.rds")
#
# ### create object to use for prediction
posterior_samples_T1D_obj <- list(post = posterior_samples_T1D$samples)
class(posterior_samples_T1D_obj) <- "T1D"
predictions_dataset_type1_with_T <- predict(posterior_samples_T1D_obj, dataset_type1, rcs_parms)

predictions_dataset_type1_with_T_new <- predict(posterior_samples_T1D_obj, dataset_type1, rcs_parms) %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, prob = 0.025), UCI = quantile(x, prob = 0.975))
  }) %>%
  bind_rows()

#MERGE TO MY_T2D
dataset_type1 <- cbind(dataset_type1, predictions_dataset_type1_with_T_new)



##T2D/NOT-EARLY-INSULIN-TREATED MODEL -----------------------------------------------------------------
posterior_samples_T2D <- readRDS("type_2_model_posteriors_thin_100.rds")
# 
posterior_samples_T2D_obj <- list(post = posterior_samples_T2D$samples)
class(posterior_samples_T2D_obj) <- "T2D"



predictions_dataset.type2_new <- predict(posterior_samples_T2D_obj, dataset_type2) %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, prob = 0.025), UCI = quantile(x, prob = 0.975))
  }) %>%
  bind_rows()


predictions_dataset_type2 <- predict(posterior_samples_T2D_obj, dataset_type2)

#MERGE TO MY_T2D
dataset_type2 <- cbind(dataset_type2, predictions_dataset.type2_new)


#merge together --------------------------------------------------------------------------
External_joint <- full_join(dataset_type1, dataset_type2)

#Plotting ------------------------------------------------------------------------------
setwd("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper")

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
auc_T1D_with_T <- calc_auroc(dataset_type1$M, predictions_dataset_type1_with_T, thinning = 10)
quantile(auc_T1D_with_T, probs = c(0.025, 0.5, 0.975)) # thinning = 10
# 2.5%       50%     97.5% 
# 0.9719786 0.9809626 0.9837487 

## Type 2 UNITED
auc_T2D_new <- calc_auroc(dataset_type2$M, predictions_dataset_type2, thinning = 10)
quantile(auc_T2D_new, probs = c(0.025, 0.5, 0.975)) # thinning = 10
# 2.5%       50%     97.5% 
# 0.8990476 0.9138095 0.9276190 

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
### Biomarker models
roc_T1D_with_T <- calc_roc(dataset_type1$M, predictions_dataset_type1_with_T, thinning = 1000)

# plot for ROC with grey being iterations, black being the ROC for average prediction
plot_roc_T1D_with_T <- ggplot() +
  ## all iterations
  geom_path(
    data = roc_T1D_with_T,
    aes(x = 1-specificities, y= sensitivities, group = iteration), colour = "grey"
  ) +
  ## average predictions
  geom_path(
    data = pROC::roc(response = dataset_type1$M, predictor = colMeans(predictions_dataset_type1_with_T)) %>%
      magrittr::extract(c(2:3)) %>%
      as.data.frame(),
    aes(x = 1-specificities, y= sensitivities), colour = "black"
  )

## Type 2 UNITED
roc_T2D <- calc_roc(dataset_type2$M, predictions_dataset_type2, thinning = 1000)

# plot for ROC with grey being iterations, black being the ROC for average prediction
plot_roc_T2D <- ggplot() +
  ## all iterations
  geom_path(
    data = roc_T2D,
    aes(x = 1-specificities, y= sensitivities, group = iteration), colour = "grey"
  ) +
  ## average predictions
  geom_path(
    data = pROC::roc(response = dataset_type2$M, predictor = colMeans(predictions_dataset_type2)) %>%
      magrittr::extract(c(2:3)) %>%
      as.data.frame(),
    aes(x = 1-specificities, y= sensitivities), colour = "black"
  )



#:------------------------------------------------------------
roc_curves <- data.frame(prob = colMeans(predictions_dataset_type1_with_T)) %>%
  cbind(Mody = dataset_type1$M) %>%
  pROC::roc(response = Mody, predictor = prob) %>%
  magrittr::extract(2:3) %>%
  as.data.frame() %>%
  mutate(
    auc =  unname(data.frame(prob = colMeans(predictions_dataset_type1_with_T)) %>%
                    cbind(Mody = dataset_type1$M) %>%
                    pROC::roc(response = Mody, predictor = prob) %>%
                    magrittr::extract(c(9)) %>%
                    unlist()),
    auc_low = quantile(auc_T1D_with_T, probs = c(0.025)),
    auc_high = quantile(auc_T1D_with_T, probs= c(0.975)),
    mean = mean(colMeans(predictions_dataset_type1_with_T), na.rm = TRUE)
  ) %>%
  mutate(Dataset = "External validation", Model = "Type 1", Calculator = "Biomarkers") %>%
  rbind( 
    data.frame(prob = colMeans(predictions_dataset_type2)) %>%
      cbind(Mody = dataset_type2$M) %>%
      pROC::roc(response = Mody, predictor = prob) %>%
      magrittr::extract(2:3) %>%
      as.data.frame() %>%
      mutate(
        auc = unname(data.frame(prob = colMeans(predictions_dataset_type2)) %>%
                       cbind(Mody = dataset_type2$M) %>%
                       pROC::roc(response = Mody, predictor = prob) %>%
                       magrittr::extract(c(9)) %>%
                       unlist()),
        auc_low = quantile(auc_T2D_new, probs = c(0.025)),
        auc_high = quantile(auc_T2D_new, probs= c(0.975)),
        mean = mean(colMeans(predictions_dataset_type2), na.rm = TRUE)
      ) %>%
      mutate(Dataset = "External validation", Model = "Type 2", Calculator = "No Biomarkers")
  )




dat_text <- roc_curves %>%
  select(-sensitivities, -specificities) %>%
  distinct() %>%
  mutate(
    auc_full = paste0("AUC: ", signif(auc, 2), " [", signif(auc_low, 2), "-", signif(auc_high, 2), "]"),
    mean = paste0("Mean prob:", signif(mean, 2)*100, "%"),
    Calculator = factor(Calculator, levels = c("Biomarkers"), labels = c("Clinical features and biomarkers"))
  )



plot_prob_external_fig3 <- patchwork::wrap_plots(
  patchwork::wrap_plots(
  # ROC T1D
  patchwork::free(
    
    roc_curves %>%
      filter(Dataset == "External validation" & Model == "Type 1") %>%
      mutate(
        Calculator = factor(Calculator, 
                            levels = c("No Biomarkers", "Biomarkers"), 
                            labels = c("Clinical features", "Clinical features and biomarkers")),
        iteration = 0
      ) %>%
      ggplot(aes(x = 1- specificities, y = sensitivities)) +
      geom_path() +
      theme_bw() +
      #facet_grid(~factor(Model, 
                         #levels = c("Type 1", "Type 2"), 
                         #labels = c("Early-insulin-treated", "Not-early-insulin-treated")), 
                 #scales = "free",
                 #labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
      scale_y_continuous("Sensitivity", labels = scales::percent) +
      scale_x_continuous("1- Specificity", labels = scales::percent) +
      theme_bw() +
      geom_label(
        data = dat_text %>%
          filter(Dataset == "External validation" & Model == "Type 1"),
        mapping = aes(x = 0.55, y = 0.1, label = auc_full, hjust = "center"),
        size = 7,
        label.r = unit(0, "pt"),
        label.padding=unit(0.4, "lines")
      ) +
      theme(
        panel.spacing.x = unit(1.5, "lines")
      )
  ),
  # ROC T2D
  patchwork::free(
    
    roc_curves %>%
      filter(Dataset == "External validation" & Model == "Type 2") %>%
      mutate(
        Calculator = factor(Calculator, 
                            levels = c("No Biomarkers", "Biomarkers"), 
                            labels = c("Clinical features", "Clinical features and biomarkers")),
        iteration = 0
      ) %>%
      ggplot(aes(x = 1- specificities, y = sensitivities)) +
      geom_path() +
      theme_bw() +
      #facet_grid(~factor(Model, 
                         #levels = c("Type 1", "Type 2"), 
                         #labels = c("Early-insulin-treated", "Not-early-insulin-treated")), 
                 #scales = "free",
                 #labeller = label_wrap_gen(width = 46, multi_line =  TRUE)) +
      scale_y_continuous("Sensitivity", labels = scales::percent) +
      scale_x_continuous("1- Specificity", labels = scales::percent) +
      theme_bw() +
      geom_label(
        data = dat_text %>%
          filter(Dataset == "External validation" & Model == "Type 2"),
        mapping = aes(x = 0.55, y = 0.1, label = auc_full, hjust = "center"),
        size = 7,
        label.r = unit(0, "pt"),
        label.padding=unit(0.4, "lines")
      ) +
      theme(
        panel.spacing.x = unit(1.5, "lines")
      )
  ),
nrow = 1, ncol = 2
  ),
  patchwork::wrap_plots(
    
    patchwork::wrap_plots(
      
      # Density
      External_joint %>%
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
        ylab("Non-MODY \n (n=1005)"),
      #point
      External_joint %>%
        filter(M == 1) %>%
        select(M, prob) %>%
        rename("Mody" = "M") %>%
        mutate(
          Mody = factor(Mody, levels = c(0, 1), labels = c("Non-MODY", "MODY")),
        ) %>%
        ggplot() +
        geom_point(aes(x = prob, y=0), position = position_jitter(height = 0.13, width = 0.1, seed = 24)) +
        geom_vline(xintercept = 0.05) +
        coord_cartesian(xlim =c(0, 1), ylim = c(-0.15, 0.15)) +
        scale_x_continuous(labels = scales::percent) +
        theme_classic() +
        theme(
          axis.ticks.y = element_blank(), 
          axis.text.y = element_blank(), 
          #axis.title.y = element_blank()
        ) +
        ylab("MODY \n cases \n (n=20)") +
        xlab("Model probabilities"),
      ncol = 1, nrow = 2, heights = c(3.5,1.5)
    )
    
  ),
  
  nrow = 2, ncol = 1
  
  
) + patchwork::plot_annotation(tag_levels = list(c("A.1", "A.2", "B", ""))) +
  patchwork::plot_layout(
    design = "
    AAAAA
    #BBB#
    "
  ) &
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 11), 
    plot.margin = margin(10,10,10,10)
  ) 


pdf("figures/Figure3.pdf", width = 10, height = 8)
plot_prob_external_fig3
dev.off()

