#:--------------------------------------------------------
# Supplementary Figure 9 
#
# In this file we plot the calibration plots for UNITED models
#
#:--------------------------------------------------------

# load libraries
library(tidyverse)
library(ggplot2)
library(patchwork)

# load files required
predictions_dataset.UNITED_type1_with_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_all_genes_with_T.rds")
predictions_dataset.UNITED_type2_new <- readRDS("model_predictions/predictions_dataset.UNITED_type2_all_genes_new.rds")


# predictions_dataset.UNITED_type1_with_T_full <- readRDS("model_predictions/predictions_dataset.UNITED_type1_all_genes_with_T_full.rds")
# predictions_dataset.UNITED_type2_new_full <- readRDS("model_predictions/predictions_dataset.UNITED_type2_all_genes_new_full.rds")


# load functions needed
source("data/create_data.R")



## Load population representative dataset
dataset.UNITED_type1 <- create_data(dataset = "united t1d", commonmody = FALSE)
# we don't do this to make the plotting easier since these patients don't have the correct biomarker status for testing (cpeptide negative or antibody positive)
# and hence their probabilities are all the same <0.001
# ## if MODY testing missing, change to 0
# mutate(M = ifelse(is.na(M), 0, M))    
dataset.UNITED_type2 <- create_data(dataset = "united t2d", commonmody = FALSE) %>%
  
  ## if MODY testing missing, change to 0
  mutate(M = ifelse(is.na(M), 0, M))

#:--------------------------------------------

## Calibration plots for UNITED type 1


#Define grouping values
#this takes the probabilities in the "probs" column and filters them based on whether their corresponding line in the dataset has a non-missing M
#i.e. only looking at the probabilities for those individuals that have the correct biomarker status (cpeptide positive AND antibody negative)
#that were therefore eligible for testing and will have a prob that isn't (and is greater than) <0.001
grouping_values_UNITED_type1_with_T <- predictions_dataset.UNITED_type1_with_T$prob[which(!is.na(dataset.UNITED_type1$M))] %>% 
  #we then find the 60th and 80th quantile for these individuals
  quantile(probs = c(0.33, 0.66))

# interim dataset combining probabilities, MODY status and grouping variable
interim_dataset <- data.frame(
  prob = predictions_dataset.UNITED_type1_with_T$prob[which(!is.na(dataset.UNITED_type1$M))],
  M = dataset.UNITED_type1$M[which(!is.na(dataset.UNITED_type1$M))]
) %>%
  mutate(
    group = as.numeric(cut(prob, breaks = c(0, grouping_values_UNITED_type1_with_T, 1)))
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
      mean = mean(predictions_dataset.UNITED_type1_with_T$prob[which(is.na(dataset.UNITED_type1$M))]),
      fit = 0, upr = NA, lwr = NA, shape = "not mody test"
    )
  )

plot_calibration_UNITED_type1_with_T <- dataset_plot %>%
  ggplot(aes(x = mean, y = fit)) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed") +
  geom_smooth(method = 'lm', formula = y ~ x,
              aes(fill = after_scale(color)), alpha = 0) +
  geom_point(aes(shape = shape), size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr)) +
  scale_shape_manual(values = c("mody tested" = "triangle", "not mody test" = "square")) +
  xlim(0, 1) +
  ylim(0, 1) +
  coord_cartesian(ylim =c(0, 0.4), xlim =c(0, 0.4)) +
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
grouping_values_UNITED_type2_new <- predictions_dataset.UNITED_type2_new$prob[which(!is.na(dataset.UNITED_type2$M))] %>% 
  #we then find the 20th and 80th quantile for these individuals
  quantile(probs = c(0.33, 0.66))


interim_dataset <- data.frame(
  prob = predictions_dataset.UNITED_type2_new$prob,
  M = dataset.UNITED_type2$M
) %>%
  mutate(
    group = as.numeric(cut(prob, breaks = c(0, grouping_values_UNITED_type2_new, 1)))
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


plot_calibration_UNITED_type2_new <- dataset_plot %>%
  ggplot(aes(x = mean, y = fit)) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed") +
  geom_smooth(method = 'lm', formula = y ~ x,
              aes(fill = after_scale(color)), alpha = 0) +
  geom_point(shape = "triangle", size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr)) +
  xlim(0, 1) +
  ylim(0, 1) +
  coord_cartesian(ylim =c(0, 0.6), xlim =c(0, 0.6)) +
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
  
  plot_calibration_UNITED_type1_with_T,
  
  plot_calibration_UNITED_type2_new,
  
  ncol = 1
  
) +
  patchwork::plot_annotation(tag_levels = 'A') +
  plot_layout(guides = 'collect')




pdf("Supplementary Material/Outputs/supfig9.pdf", width = 6, height = 9)
plot_calibration
dev.off()
