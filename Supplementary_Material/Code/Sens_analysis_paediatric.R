######################################################################################
#Sensitivity Analysis: Paediatric UNITED

#######################################################################################
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
# load datasets
## Load population representative dataset
#readin MYDIABETES
#load T1D
load("Data/MY_T1D.RData")


# MyDiabetes ------------------------------------------------------------------------------------
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


# UNITED paediatric ---------------------------------------------------------------------
dataset.UNITED_p <- create_data(dataset = "united t1d pediatrics", 
                                commonmody = FALSE, 
                                id = TRUE)

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
  select(MY_ID, study, agerec, agedx, sex, bmi, pardm, 
         insoroha, hba1c, C, A, M, Gene, biomark_status) 

dataset.UNITED_type1p <- dataset.UNITED_type1p %>%
  mutate(study = "UNITED paediatric") %>%
  filter(!is.na(T))
dataset_type1 <- full_join(dataset.MYDIABETES_type1, 
                           dataset.UNITED_type1p, 
                           by = c("study","agerec", "agedx", "sex", "bmi", 
                                  "pardm", "insoroha", "hba1c", "C", "A", 
                                  "M", "Gene")) %>%
  mutate(M = ifelse(is.na(M), 0, M), 
         Gene = ifelse(Gene == "", NA, Gene), 
         biomark_status = ifelse(is.na(biomark_status), 
                                 ifelse(C == 1 & A == 0, 0, 1),
                                 biomark_status))
dataset_type2 <- dataset_type1 %>%
  rename(agerec_original = agerec) %>%
  mutate(agerec = ifelse(agerec_original > 50, 50, agerec_original)) %>%
  filter(agerec < 18)
dataset_type1 <- dataset_type1 %>%
  rename(agerec_original = agerec) %>%
  mutate(agerec = ifelse(agerec_original > 50, 50, agerec_original)) %>%
  filter(agedx < 18)


table(dataset_type1$M)
table(dataset_type2$M)
#run predictions ------------------------------------------------------------------------------
## In MYDIABETES
setwd("New_Data_Predictions")
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

predictions_dataset_type1_with_T2 <- predict(posterior_samples_T1D_obj, dataset_type2, rcs_parms)

predictions_dataset_type1_with_T_new2 <- predict(posterior_samples_T1D_obj, dataset_type2, rcs_parms) %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, prob = 0.025), UCI = quantile(x, prob = 0.975))
  }) %>%
  bind_rows()

#MERGE TO MY_T2D
dataset_type1 <- cbind(dataset_type1, predictions_dataset_type1_with_T_new)
dataset_type2 <- cbind(dataset_type2, predictions_dataset_type1_with_T_new2)




#Figure prep --------------------------------------------------------------------------------------

setwd("..")

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

auc_T1D_with_T2 <- calc_auroc(dataset_type2$M, predictions_dataset_type1_with_T2, thinning = 10)
quantile(auc_T1D_with_T2, probs = c(0.025, 0.5, 0.975)) # thinning = 10



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
roc_T1D_with_T <- calc_roc(dataset_type1$M, 
                           predictions_dataset_type1_with_T, 
                           thinning = 1000)
roc_T1D_with_T2 <- calc_roc(dataset_type2$M, 
                           predictions_dataset_type1_with_T2, 
                           thinning = 1000)

# plot for ROC with grey being iterations, black being the ROC for average prediction
plot_roc_T1D_with_T <- ggplot() +
  ## all iterations
  geom_path(
    data = roc_T1D_with_T,
    aes(x = 1-specificities, y= sensitivities, group = iteration), colour = "grey"
  ) +
  ## average predictions
  geom_path(
    data = pROC::roc(response = dataset_type1$M, 
                     predictor = colMeans(predictions_dataset_type1_with_T)) %>%
      magrittr::extract(c(2:3)) %>%
      as.data.frame(),
    aes(x = 1-specificities, y= sensitivities), colour = "black"
  )
plot_roc_T1D_with_T2 <- ggplot() +
  ## all iterations
  geom_path(
    data = roc_T1D_with_T2,
    aes(x = 1-specificities, y= sensitivities, group = iteration), colour = "grey"
  ) +
  ## average predictions
  geom_path(
    data = pROC::roc(response = dataset_type2$M, 
                     predictor = colMeans(predictions_dataset_type1_with_T2)) %>%
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
  mutate(Dataset = "External validation", Model = "Type 1", Calculator = "Biomarkers") 
roc_curves2 <- data.frame(prob = colMeans(predictions_dataset_type1_with_T2)) %>%
  cbind(Mody = dataset_type2$M) %>%
  pROC::roc(response = Mody, predictor = prob) %>%
  magrittr::extract(2:3) %>%
  as.data.frame() %>%
  mutate(
    auc =  unname(data.frame(prob = colMeans(predictions_dataset_type1_with_T2)) %>%
                    cbind(Mody = dataset_type2$M) %>%
                    pROC::roc(response = Mody, predictor = prob) %>%
                    magrittr::extract(c(9)) %>%
                    unlist()),
    auc_low = quantile(auc_T1D_with_T2, probs = c(0.025)),
    auc_high = quantile(auc_T1D_with_T2, probs= c(0.975)),
    mean = mean(colMeans(predictions_dataset_type1_with_T2), na.rm = TRUE)
  ) %>%
  mutate(Dataset = "External validation", Model = "Type 1", Calculator = "Biomarkers") 




dat_text <- roc_curves %>%
  select(-sensitivities, -specificities) %>%
  distinct() %>%
  mutate(
    auc_full = paste0("AUC: ", signif(auc, 2), " [", signif(auc_low, 2), "-", signif(auc_high, 2), "]"),
    mean = paste0("Mean prob:", signif(mean, 2)*100, "%"),
    Calculator = factor(Calculator, levels = c("Biomarkers"), labels = c("Clinical features and biomarkers"))
  )
dat_text2 <- roc_curves2 %>%
  select(-sensitivities, -specificities) %>%
  distinct() %>%
  mutate(
    auc_full = paste0("AUC: ", signif(auc, 2), " [", signif(auc_low, 2), "-", signif(auc_high, 2), "]"),
    mean = paste0("Mean prob:", signif(mean, 2)*100, "%"),
    Calculator = factor(Calculator, levels = c("Biomarkers"), labels = c("Clinical features and biomarkers"))
  )





# ## Calibration plots for UNITED type 1
# 
# 
# #Define grouping values
# #this takes the probabilities in the "probs" column and filters them based on whether their corresponding line in the dataset has a non-missing M
# #i.e. only looking at the probabilities for those individuals that have the correct biomarker status (cpeptide positive AND antibody negative)
# #that were therefore eligible for testing and will have a prob that isn't (and is greater than) <0.001
# grouping_values_type1_with_T <- predictions_dataset_type1_with_T_new$prob[which(!is.na(dataset_type1$M))] %>% 
#   #we then find the 60th and 80th quantile for these individuals
#   quantile(probs = c(0.33, 0.66))
# 
# # interim dataset combining probabilities, MODY status and grouping variable
# interim_dataset <- data.frame(
#   prob = predictions_dataset_type1_with_T_new$prob[which(!is.na(dataset_type1$M))],
#   M = dataset_type1$M[which(!is.na(dataset_type1$M))]
# ) %>%
#   mutate(
#     group = as.numeric(cut(prob, breaks = c(0, grouping_values_type1_with_T, 1)))
#   )
# 
# # dataset with plot information
# dataset_plot <- NULL
# 
# for (i in 1:length(unique(interim_dataset$group))) {
#   
#   # select entries for this group
#   interim_data <- interim_dataset %>%
#     filter(group == i) %>%
#     mutate(M = factor(M))
#   
#   # fit the model  calculating the CI for each point
#   interim_model <- glm(M ~ 1, data = interim_data, 
#                        family=binomial(link='logit'))
#   
#   # calculate the confidence interval
#   prediction <- predict(interim_model, 
#                         newdata = data.frame(prob = mean(interim_data$prob)), 
#                         type = "link", 
#                         se.fit = TRUE)
#   critval <- 1.96 ## approx 95% CI
#   upr <- prediction$fit + (critval * prediction$se.fit)
#   lwr <- prediction$fit - (critval * prediction$se.fit)
#   fit <- prediction$fit
#   
#   # turn values into a probability prediction (since above is the link function value)
#   fit2 <- interim_model$family$linkinv(fit)
#   upr2 <- interim_model$family$linkinv(upr)
#   lwr2 <- interim_model$family$linkinv(lwr)
#   
#   # combine into plot information dataset
#   dataset_plot <- rbind(
#     dataset_plot,
#     data.frame(mean = mean(interim_data$prob), 
#                fit = fit2, 
#                upr = upr2, 
#                lwr = lwr2)
#   )
#   
# }
# 
# # combine with value for C-/A+ patients
# dataset_plot <- dataset_plot %>%
#   mutate(shape = "mody tested") %>%
#   rbind(
#     data.frame(
#       mean = mean(predictions_dataset_type1_with_T_new$prob[which(is.na(dataset_type1$M))]),
#       fit = 0, upr = NA, lwr = NA, shape = "not mody test"
#     )
#   )
# 
# plot_calibration_type1_with_T <- dataset_plot %>%
#   ggplot(aes(x = mean, y = fit)) +
#   geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed") +
#   geom_smooth(method = 'lm', formula = y ~ x,
#               aes(fill = after_scale(color)), alpha = 0) +
#   geom_point(aes(shape = shape), size = 3) +
#   geom_errorbar(aes(ymin = lwr, ymax = upr)) +
#   scale_shape_manual(values = c("mody tested" = "triangle", 
#                                 "not mody test" = "square")) +
#   xlim(0, 1) +
#   ylim(0, 1) +
#   coord_cartesian(ylim =c(0, 1), xlim =c(0, 1)) +
#   scale_y_continuous(labels = scales::percent) +
#   scale_x_continuous(labels = scales::percent) +
#   xlab("Model predictions") +
#   ylab("Observed probability") +
#   theme_light() +
#   theme(
#     axis.text = element_text(size = 12),
#     strip.text = element_text(size = 16),
#     title = element_text(size = 17),
#     legend.position = "none",
#     legend.text = element_text(size = 12))




plot_prob_SENS_PAED <- patchwork::wrap_plots(
    # ROC T1D
  roc_curves %>%
    filter(Dataset == "External validation" & Model == "Type 1") %>%
    mutate(
      Calculator = factor(Calculator, 
                          levels = c("Biomarkers"), 
                          labels = c("Clinical features and biomarkers")),
      iteration = 0
      ) %>%
    ggplot(aes(x = 1- specificities, y = sensitivities)) +
    geom_path() +
    theme_bw() +
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
      ),
  patchwork::wrap_plots(
    # Density
    dataset_type1 %>%
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
      ylab("Non-MODY \n (n=715)"),
    #point
    dataset_type1 %>%
      filter(M == 1) %>%
      select(M, prob) %>%
      rename("Mody" = "M") %>%
      mutate(
        Mody = factor(Mody, levels = c(0, 1), labels = c("Non-MODY", "MODY")),
        ) %>%
      ggplot() +
      geom_point(aes(x = prob, y=0), 
                 #position = position_jitter(height = 0.13, width = 0.1, seed = 24)
                 ) +
      geom_vline(xintercept = 0.05) +
      coord_cartesian(xlim =c(0, 1), ylim = c(-0.15, 0.15)) +
      scale_x_continuous(labels = scales::percent) +
      theme_classic() +
      theme(
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        #axis.title.y = element_blank()
        ) +
      ylab("MODY \n cases \n (n=3)") +
      xlab("Model probabilities"),
  ncol = 1, nrow = 2, heights = c(3.5,1.5)
  ),
  roc_curves2 %>%
    filter(Dataset == "External validation" & Model == "Type 1") %>%
    mutate(
      Calculator = factor(Calculator, 
                          levels = c("Biomarkers"), 
                          labels = c("Clinical features and biomarkers")),
      iteration = 0
    ) %>%
    ggplot(aes(x = 1- specificities, y = sensitivities)) +
    geom_path() +
    theme_bw() +
    scale_y_continuous("Sensitivity", labels = scales::percent) +
    scale_x_continuous("1- Specificity", labels = scales::percent) +
    theme_bw() +
    geom_label(
      data = dat_text2 %>%
        filter(Dataset == "External validation" & Model == "Type 1"),
      mapping = aes(x = 0.55, y = 0.1, label = auc_full, hjust = "center"),
      size = 7,
      label.r = unit(0, "pt"),
      label.padding=unit(0.4, "lines")
    ) +
    theme(
      panel.spacing.x = unit(1.5, "lines")
    ),
  patchwork::wrap_plots(
    # Density
    dataset_type2 %>%
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
      ylab("Non-MODY \n (n=436)"),
    #point
    dataset_type2 %>%
      filter(M == 1) %>%
      select(M, prob) %>%
      rename("Mody" = "M") %>%
      mutate(
        Mody = factor(Mody, levels = c(0, 1), labels = c("Non-MODY", "MODY")),
      ) %>%
      ggplot() +
      geom_point(aes(x = prob, y=0), 
                 #position = position_jitter(height = 0.13, #width = 0.1, seed = 25)
                 ) +
      geom_vline(xintercept = 0.05) +
      coord_cartesian(xlim =c(0, 1), ylim = c(-0.15, 0.15)) +
      scale_x_continuous(labels = scales::percent) +
      theme_classic() +
      theme(
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        #axis.title.y = element_blank()
      ) +
      ylab("MODY \n cases \n (n=2)") +
      xlab("Model probabilities"),
    ncol = 1, nrow = 2, heights = c(3.5,1.5)
  ),
nrow = 2, ncol = 2
) + patchwork::plot_annotation(tag_levels = list(c("A", "B", " ","C", "D")))  &
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 11), 
    plot.margin = margin(10,10,10,10)
  ) 


pdf("Figures/Supp_PAED.pdf", width = 20, height = 18)
plot_prob_SENS_PAED
dev.off()

dataset_type1 %>%
  filter(M==1) %>%
  select(Gene, prob, agedx, agerec, C, A)
table(dataset_type1$M)

dataset_type2 %>%
  filter(M==1) %>%
  select(Gene, prob, agedx, agerec, C, A)
table(dataset_type2$M)


PAED_5PERC <- dataset_type1 %>%
  summarise(totalover = sum(prob >= 0.05),
            ncasespickedup = sum(prob >= 0.05 & M ==1),
            PPV = (sum(prob >= 0.05 & M ==1)/sum(prob >= 0.05))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.05 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.05 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.05 & M == 0)/sum(prob < 0.05))*100,
            Sensitivity = (sum(prob >= 0.05 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.05 & M == 0)/sum(M==0))*100)

PAED_5PERC2 <- dataset_type2 %>%
  summarise(totalover = sum(prob >= 0.05),
            ncasespickedup = sum(prob >= 0.05 & M ==1),
            PPV = (sum(prob >= 0.05 & M ==1)/sum(prob >= 0.05))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.05 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.05 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.05 & M == 0)/sum(prob < 0.05))*100,
            Sensitivity = (sum(prob >= 0.05 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.05 & M == 0)/sum(M==0))*100)
