#:--------------------------------------------------------
# Supplementary Figure 15
#
# In this file we plot the calibration plots for MODY models in the External validation
#
#:--------------------------------------------------------
setwd("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper")

# load libraries -------------------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(patchwork)
library(readxl)
library(writexl)


# load functions needed --------------------------------------------------------------------
source("new_data_predictions/prediction_functions.R")
source("data/create_data.R")

#readin datasets --------------------------------------------------------------------------------------
# MyDiabetes ------------------------------------------------------------------------------------
#load T1D
load("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper/MY_T1D.RData")
#load T2D 
load("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper/MY_T2D.RData")
# LIMIT TO ONLY WHITE ETHNICITY
MY_T1D <- MY_T1D %>%
  filter(ethnicity == "White") %>%
  mutate(M = ifelse(biomark_status == 0, M, NA))

table(MY_T1D$M)
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
  mutate(
         Gene = ifelse(Gene == "", NA, Gene), 
         biomark_status = ifelse(is.na(biomark_status), 
                                 ifelse(C == 1 & A == 0, 0, 1),
                                 biomark_status)) 
#dataset_type1 <- dataset_type1 %>%
  #rename(agerec_original = agerec) %>%
  #mutate(agerec = ifelse(agerec_original > 50, 50, agerec_original))

# Not-early-insulin-treated
dataset.MYDIABETES_type2 <- MY_T2D %>%
  mutate(study = "MYDIABETES") %>%
  select(MY_ID, study, agerec, agedx, sex, bmi, pardm, insoroha, hba1c, C, A, M, Gene, biomark_status) 

dataset.UNITED_type2p <- dataset.UNITED_type2p %>%
  mutate(study = "UNITED paediatric") %>%
  select(id, study, agerec, agedx, sex, bmi, pardm, insoroha, hba1c, C, A, M, Gene) 

dataset_type2 <- full_join(dataset.MYDIABETES_type2, dataset.UNITED_type2p, by = c("study","agerec", "agedx", "sex", "bmi", "pardm", "insoroha", "hba1c", "C", "A", "M", "Gene")) %>%
  mutate(
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
predictions_dataset_type1_with_T_new <- predict(posterior_samples_T1D_obj, dataset_type1, rcs_parms)

predictions_dataset_type1_with_T <- predict(posterior_samples_T1D_obj, dataset_type1, rcs_parms) %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, prob = 0.025), UCI = quantile(x, prob = 0.975))
  }) %>%
  bind_rows()

#MERGE TO MY_T2D
dataset_type1 <- cbind(dataset_type1, predictions_dataset_type1_with_T)



##T2D/NOT-EARLY-INSULIN-TREATED MODEL -----------------------------------------------------------------
posterior_samples_T2D <- readRDS("type_2_model_posteriors_thin_100.rds")
# 
posterior_samples_T2D_obj <- list(post = posterior_samples_T2D$samples)
class(posterior_samples_T2D_obj) <- "T2D"



predictions_dataset_type2_new <- predict(posterior_samples_T2D_obj, dataset_type2) %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, prob = 0.025), UCI = quantile(x, prob = 0.975))
  }) %>%
  bind_rows()


predictions_dataset_type2 <- predict(posterior_samples_T2D_obj, dataset_type2)

#MERGE TO MY_T2D
dataset_type2 <- cbind(dataset_type2, predictions_dataset_type2_new)



#:--------------------------------------------
setwd("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper")



## Calibration plots for UNITED type 1


#Define grouping values
#this takes the probabilities in the "probs" column and filters them based on whether their corresponding line in the dataset has a non-missing M
#i.e. only looking at the probabilities for those individuals that have the correct biomarker status (cpeptide positive AND antibody negative)
#that were therefore eligible for testing and will have a prob that isn't (and is greater than) <0.001
grouping_values_type1_with_T <- predictions_dataset_type1_with_T$prob[which(!is.na(dataset_type1$M))] %>% 
  #we then find the 60th and 80th quantile for these individuals
  quantile(probs = c(0.33, 0.66))

# interim dataset combining probabilities, MODY status and grouping variable
interim_dataset <- data.frame(
  prob = predictions_dataset_type1_with_T$prob[which(!is.na(dataset_type1$M))],
  M = dataset_type1$M[which(!is.na(dataset_type1$M))]
) %>%
  mutate(
    group = as.numeric(cut(prob, breaks = c(0, grouping_values_type1_with_T, 1)))
  )

# dataset with plot information
dataset_plot <- NULL

for (i in 1:length(unique(interim_dataset$group))) {
  
  # select entries for this group
  interim_data <- interim_dataset %>%
    filter(group == i) %>%
    mutate(M = factor(M))
  
  # fit the model  calculating the CI for each point
  interim_model <- glm(M ~ 1, data = interim_data, family=binomial(link='logit'))
  
  # calculate the confidence interval
  prediction <- predict(interim_model, newdata = data.frame(prob = mean(interim_data$prob)), type = "link", se.fit = TRUE)
  critval <- 1.96 ## approx 95% CI
  upr <- prediction$fit + (critval * prediction$se.fit)
  lwr <- prediction$fit - (critval * prediction$se.fit)
  fit <- prediction$fit
  
  # turn values into a probability prediction (since above is the link function value)
  fit2 <- interim_model$family$linkinv(fit)
  upr2 <- interim_model$family$linkinv(upr)
  lwr2 <- interim_model$family$linkinv(lwr)
  
  # combine into plot information dataset
  dataset_plot <- rbind(
    dataset_plot,
    data.frame(mean = mean(interim_data$prob), fit = fit2, upr = upr2, lwr = lwr2)
  )
  
}

# combine with value for C-/A+ patients
dataset_plot <- dataset_plot %>%
  mutate(shape = "mody tested") %>%
  rbind(
    data.frame(
      mean = mean(predictions_dataset_type1_with_T$prob[which(is.na(dataset_type1$M))]),
      fit = 0, upr = NA, lwr = NA, shape = "not mody test"
    )
  )

plot_calibration_type1_with_T <- dataset_plot %>%
  ggplot(aes(x = mean, y = fit)) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed") +
  geom_smooth(method = 'lm', formula = y ~ x,
              aes(fill = after_scale(color)), alpha = 0) +
  geom_point(aes(shape = shape), size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr)) +
  scale_shape_manual(values = c("mody tested" = "triangle", "not mody test" = "square")) +
  xlim(0, 1) +
  ylim(0, 1) +
  coord_cartesian(ylim =c(0, 1), xlim =c(0, 1)) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  xlab("Model predictions") +
  ylab("Observed probability") +
  theme_light() +
  theme(
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 16),
    title = element_text(size = 17),
    legend.position = "none",
    legend.text = element_text(size = 12))



#:--------------------------------------------

## Calibration plots for UNITED type 2


#Define grouping values
#this takes the probabilities in the "probs" column and filters them based on whether their corresponding line in the dataset has a non-missing M
#i.e. only looking at the probabilities for those individuals weren't tested for MODY
#in this case this is irrelevant as there are none and we altered any missing to 0 
#therefore looking for quantiles on all predicted probs
grouping_values_type2_new <- predictions_dataset_type2_new$prob[which(!is.na(dataset_type2$M))] %>% 
  #we then find the 20th and 80th quantile for these individuals
  quantile(probs = c(0.33, 0.66))


interim_dataset <- data.frame(
  prob = predictions_dataset_type2_new$prob,
  M = dataset_type2$M
) %>%
  mutate(
    group = as.numeric(cut(prob, breaks = c(0, grouping_values_type2_new, 1)))
  )

dataset_plot <- NULL

for (i in 1:length(unique(interim_dataset$group))) {
  
  # select entries for this group
  interim_data <- interim_dataset %>%
    filter(group == i) %>%
    mutate(M = factor(M))
  
  # fit the model  calculating the CI for each point
  interim_model <- glm(M ~ 1, data = interim_data, family=binomial(link='logit'))
  
  # calculate the confidence interval
  prediction <- predict(interim_model, newdata = data.frame(prob = mean(interim_data$prob)), type = "link", se.fit = TRUE)
  critval <- 1.96 ## approx 95% CI
  upr <- prediction$fit + (critval * prediction$se.fit)
  lwr <- prediction$fit - (critval * prediction$se.fit)
  fit <- prediction$fit
  
  # turn values into a probability prediction (since above is the link function value)
  fit2 <- interim_model$family$linkinv(fit)
  upr2 <- interim_model$family$linkinv(upr)
  lwr2 <- interim_model$family$linkinv(lwr)
  
  # combine into plot information dataset
  dataset_plot <- rbind(
    dataset_plot,
    data.frame(mean = mean(interim_data$prob), fit = fit2, upr = upr2, lwr = lwr2)
  )
  
}


plot_calibration_external_val_type2_new <- dataset_plot %>%
  ggplot(aes(x = mean, y = fit)) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed") +
  geom_smooth(method = 'lm', formula = y ~ x,
              aes(fill = after_scale(color)), alpha = 0) +
  geom_point(shape = "triangle", size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr)) +
  xlim(0, 1) +
  ylim(0, 1) +
  coord_cartesian(ylim =c(0, 1), xlim =c(0, 1)) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  xlab("Model predictions") +
  ylab("Observed probability") +
  theme_light() +
  theme(
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 16),
    title = element_text(size = 17),
    legend.position = "none",
    legend.text = element_text(size = 12))





plot_calibration <- patchwork::wrap_plots(
  
  plot_calibration_type1_with_T,
  
  plot_calibration_external_val_type2_new,
  
  ncol = 1
  
) +
  patchwork::plot_annotation(tag_levels = 'A') +
  plot_layout(guides = 'collect') &
  theme(
    legend.position = "bottom"
  )

pdf("Supplementary Material/Outputs/supfig15.pdf", width = 6, height = 9)
plot_calibration
dev.off()



