#:--------------------------------------------------------
#
# In this file we plot the calibration plots for UNITED models
#
#:--------------------------------------------------------

# load libraries
library(tidyverse)
library(ggplot2)
library(patchwork)

# load files required
predictions_dataset.UNITED_type1_all_genes_with_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_all_genes_with_T.rds")
predictions_dataset.UNITED_type2_all_genes_new <- readRDS("model_predictions/predictions_dataset.UNITED_type2_all_genes_new.rds")


predictions_dataset.UNITED_type1_all_genes_with_T_full <- readRDS("model_predictions/predictions_dataset.UNITED_type1_all_genes_with_T_full.rds")
predictions_dataset.UNITED_type2_all_genes_new_full <- readRDS("model_predictions/predictions_dataset.UNITED_type2_all_genes_new_full.rds")


# load functions needed
source("data/create_data.R")



## Load population representative dataset
dataset.UNITED_type1_all_genes <- create_data(dataset = "united t1d", commonmody = FALSE)
# ## if MODY testing missing, change to 0
# mutate(M = ifelse(is.na(M), 0, M))    # we don't do this to make the plotting easier since these patients don't have biomarkers and hence their probability is <0.001

dataset.UNITED_type2_all_genes <- create_data(dataset = "united t2d", commonmody = FALSE) %>%
  
  ## if MODY testing missing, change to 0
  mutate(M = ifelse(is.na(M), 0, M))



#:--------------------------------------------

## Calibration plots for UNITED type 1





legend_dataset <- cbind(
  x = c(-1, -1, -1, -1),
  y = c(-1, -1, -1, -1),
  alpha = c(0.1, 0.2, 0.3, 0.4)
) %>%
  as.data.frame()

grouping_values_UNITED_type1_all_genes_with_T <- predictions_dataset.UNITED_type1_all_genes_with_T$prob[which(!is.na(dataset.UNITED_type1_all_genes$M))] %>% quantile(probs = c(0.25, 0.5, 0.75))

# Create the points needed for the plot
## Calculate points for patients with MODY tested
## Calculate points for patients without MODY tested (assumed negative)
### Done separately because there is such a big amount of people not tested
plot_calibration_UNITED_type1_all_genes_with_T <- predictions_dataset.UNITED_type1_all_genes_with_T_full[,which(!is.na(dataset.UNITED_type1_all_genes$M))] %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  gather() %>%
  cbind(
    group = predictions_dataset.UNITED_type1_all_genes_with_T$prob[which(!is.na(dataset.UNITED_type1_all_genes$M))],
    M = dataset.UNITED_type1_all_genes$M[!is.na(dataset.UNITED_type1_all_genes$M)]
  ) %>%
  mutate(
    group = ifelse(group < grouping_values_UNITED_type1_all_genes_with_T[1], 1,
                   ifelse(group < grouping_values_UNITED_type1_all_genes_with_T[2], 2,
                          ifelse(group < grouping_values_UNITED_type1_all_genes_with_T[3], 3,
                                 ifelse(group < grouping_values_UNITED_type1_all_genes_with_T[4], 4, 5))))
  ) %>%
  group_by(group) %>%
  mutate(
    x = quantile(value, probs = c(0.5)),
    y = sum(M)/n(),
    xmin = quantile(value, probs = c(0.025)),
    xmin_75 = quantile(value, probs = c(0.125)),
    xmin_50 = quantile(value, probs = c(0.25)),
    xmin_25 = quantile(value, probs = c(0.375)),
    xmax_25 = quantile(value, probs = c(0.625)),
    xmax_50 = quantile(value, probs = c(0.75)),
    xmax_75 = quantile(value, probs = c(0.875)),
    xmax = quantile(value, probs = c(0.975))
  ) %>%
  ungroup() %>%
  select(-key, -value, -group, -M) %>%
  distinct() %>%
  rbind(
    # only plotting one patient since they are all the same
    predictions_dataset.UNITED_type1_all_genes_with_T_full[,which(is.na(dataset.UNITED_type1_all_genes$M))[1]] %>%
      as.data.frame() %>%
      gather() %>%
      cbind(
        M = dataset.UNITED_type1_all_genes$M[is.na(dataset.UNITED_type1_all_genes$M)][1]
      ) %>%
      mutate(
        M = 0 # need to do this because M is missing
      ) %>%
      mutate(
        x = quantile(value, probs = c(0.5)),
        y = sum(M)/n(),
        xmin = quantile(value, probs = c(0.025)),
        xmin_75 = quantile(value, probs = c(0.125)),
        xmin_50 = quantile(value, probs = c(0.25)),
        xmin_25 = quantile(value, probs = c(0.375)),
        xmax_25 = quantile(value, probs = c(0.625)),
        xmax_50 = quantile(value, probs = c(0.75)),
        xmax_75 = quantile(value, probs = c(0.875)),
        xmax = quantile(value, probs = c(0.975))
      ) %>%
      select(-key, -value, -M) %>%
      distinct() 
  ) %>%
  ggplot(aes(x = x, y = y)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_ribbon(aes(xmin = xmin, xmax = xmax), alpha = 0.1) +
  geom_ribbon(aes(xmin = xmin_75, xmax = xmax_75), alpha = 0.2) +
  geom_ribbon(aes(xmin = xmin_50, xmax = xmax_50), alpha = 0.3) +
  geom_ribbon(aes(xmin = xmin_25, xmax = xmax_25), alpha = 0.4) +
  geom_point(data = legend_dataset, aes(x = x, y = y, alpha = alpha), shape = "square", size = 4) +
  geom_smooth(method = 'lm', formula = y ~ x,
              aes(fill = after_scale(color)), alpha = 0) +
  geom_point(shape = "triangle") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_light() +
  xlab("Model predictions") +
  ylab("Observed probability") +
  guides(alpha = guide_legend(title = "Credible interval")) +
  scale_alpha_continuous(labels = c("95%", "75%", "50%", "25%"), range = c(0.1, 0.4)) +
  coord_cartesian(ylim =c(0, 1), xlim =c(0, 1)) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  theme(
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 16),
    title = element_text(size = 17),
    legend.position = "bottom",
    legend.text = element_text(size = 12))



grouping_values_UNITED_type2_all_genes_new <- predictions_dataset.UNITED_type2_all_genes_new$prob[which(!is.na(dataset.UNITED_type2_all_genes$M))] %>% quantile(probs = c(0.25, 0.5, 0.75))



plot_calibration_UNITED_type2_all_genes_new <- predictions_dataset.UNITED_type2_all_genes_new_full %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  gather() %>%
  cbind(
    group = predictions_dataset.UNITED_type2_all_genes_new$prob,
    M = dataset.UNITED_type2_all_genes$M
  ) %>%
  mutate(
    group = ifelse(group < grouping_values_UNITED_type2_all_genes_new[1], 1,
                   ifelse(group < grouping_values_UNITED_type2_all_genes_new[2], 2,
                          ifelse(group < grouping_values_UNITED_type2_all_genes_new[3], 3,
                                 ifelse(group < grouping_values_UNITED_type2_all_genes_new[4], 4, 5))))
  ) %>%
  mutate(group = factor(group, levels = c(1, 2, 3, 4, 5))) %>%
  group_by(group) %>%
  mutate(
    x = quantile(value, probs = c(0.5)),
    y = sum(M)/n(),
    xmin = quantile(value, probs = c(0.025)),
    xmin_75 = quantile(value, probs = c(0.125)),
    xmin_50 = quantile(value, probs = c(0.25)),
    xmin_25 = quantile(value, probs = c(0.375)),
    xmax_25 = quantile(value, probs = c(0.625)),
    xmax_50 = quantile(value, probs = c(0.75)),
    xmax_75 = quantile(value, probs = c(0.875)),
    xmax = quantile(value, probs = c(0.975))
  ) %>%
  ungroup() %>%
  select(-key, -value, -group, -M) %>%
  distinct() %>%
  ggplot(aes(x = x, y = y)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_ribbon(aes(xmin = xmin, xmax = xmax, ymin = y, ymax = y), alpha = 0.1, orientation = "y") +
  geom_ribbon(aes(xmin = xmin_75, xmax = xmax_75, ymin = y, ymax = y), alpha = 0.2, orientation = "y") +
  geom_ribbon(aes(xmin = xmin_50, xmax = xmax_50, ymin = y, ymax = y), alpha = 0.3, orientation = "y") +
  geom_ribbon(aes(xmin = xmin_25, xmax = xmax_25, ymin = y, ymax = y), alpha = 0.4, orientation = "y") +
  geom_point(data = legend_dataset, aes(x = x, y = y, alpha = alpha), shape = "square", size = 4) +
  geom_smooth(method = 'lm', formula = y ~ x,
              aes(fill = after_scale(color)), alpha = 0) +
  geom_point(shape = "triangle") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_light() +
  xlab("Model predictions") +
  ylab("Observed probability") +
  guides(alpha = guide_legend(title = "Credible interval")) +
  scale_alpha_continuous(labels = c("95%", "75%", "50%", "25%"), range = c(0.1, 0.4)) +
  coord_cartesian(ylim =c(0, 1), xlim =c(0, 1)) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  theme(
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 16),
    title = element_text(size = 17),
    legend.position = "bottom",
    legend.text = element_text(size = 12))





plot_calibration <- patchwork::wrap_plots(
  
  plot_calibration_UNITED_type1_all_genes_with_T,
  
  plot_calibration_UNITED_type2_all_genes_new,
  
  ncol = 1
  
) +
  patchwork::plot_annotation(tag_levels = 'A') +
  plot_layout(guides = 'collect') &
  theme(
    legend.position = "bottom"
  )

pdf("figures/united_calibration_plots.pdf", width = 6, height = 9)
plot_calibration
dev.off()
