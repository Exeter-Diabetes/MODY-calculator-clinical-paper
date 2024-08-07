#:--------------------------------------------------------
#   
# In this file we check Rhat and calibration plots for 
#   the riley simulation
#
#:--------------------------------------------------------

# load libraries
library(patchwork)
library(nimble)
library(tidyverse)

# load functions needed
source("data/create_data.R")
source("new_data_predictions/prediction_functions.R")

# load rds objects
bootstrap_t1d <- readRDS("bootstrap_riley_simulation/output/simulation_t1d.rds")
bootstrap_t2d <- readRDS("bootstrap_riley_simulation/output/simulation_t2d.rds")
predictions_dataset.UNITED_type1_all_genes_with_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_all_genes_with_T.rds")
predictions_dataset.UNITED_type2_all_genes_new <- readRDS("model_predictions/predictions_dataset.UNITED_type2_all_genes_new.rds")

# load datasets
## Load population representative dataset
dataset.UNITED_type1_all_genes <- create_data(dataset = "united t1d", commonmody = FALSE) %>%
  
  ## if MODY testing missing, change to 0
  mutate(M = ifelse(is.na(M), 0, M))

dataset.UNITED_type2_all_genes <- create_data(dataset = "united t2d", commonmody = FALSE)

############################################################

# Type 1 model

# Calculate Rhat values for all parameters of the model in all simulations
rhat_summary_t1d <- NULL

## append below the rhat values for each simulation
for (i in 1:bootstrap_t1d$iterations) {
  
  rhat_summary_t1d <- rbind(
    rhat_summary_t1d,
    bootstrap_t1d$simulations[[i]]$rhat_recalibration
  )
  
}

## create summary statistics
rhat_summary_t1d <- rhat_summary_t1d %>%
  
  # group by each parameters of the model
  group_by(key) %>%
  
  # create variables
  mutate(
    low_ci = quantile(rhat, probs = 0.025),
    median = quantile(rhat, probs = 0.5),
    high_ci = quantile(rhat, probs = 0.975)
  ) %>%
  
  # remove grouping
  ungroup() %>%
  
  # remove rhat variable with values
  select(-rhat) %>%
  
  # keep only unique entries
  unique()

############################

# Type 2 model

# Calculate Rhat values for all parameters of the model in all simulations
rhat_summary_t2d <- NULL

## append below the rhat values for each simulation
for (i in 1:bootstrap_t2d$iterations) {
  
  rhat_summary_t2d <- rbind(
    rhat_summary_t2d,
    bootstrap_t2d$simulations[[i]]$rhat_recalibration
  )
  
}

## create summary statistics
rhat_summary_t2d <- rhat_summary_t2d %>%
  
  # group by each parameters of the model
  group_by(key) %>%
  
  # create variables
  mutate(
    low_ci = quantile(rhat, probs = 0.025),
    median = quantile(rhat, probs = 0.5),
    high_ci = quantile(rhat, probs = 0.975)
  ) %>%
  
  # remove grouping
  ungroup() %>%
  
  # remove rhat variable with values
  select(-rhat) %>%
  
  # keep only unique entries
  unique()



############################################################

# Calibration plots for the predictions

# Type 1 model
dataset_t1d_real <- data.frame(
  id = 1:nrow(dataset.UNITED_type1_all_genes),
  recalibration_real = predictions_dataset.UNITED_type1_all_genes_with_T$prob
)

## appending the probabilities together
comparing_t1d_simulation <- NULL

for (i in 1:bootstrap_t1d$iterations) {
  
  comparing_t1d_simulation <- rbind(
    comparing_t1d_simulation,
    bootstrap_t1d$simulations[[i]]$bootstrap_dataset
  )
}

combined_probabilities_t1d <- comparing_t1d_simulation %>%
  
  # join real probability from original model
  left_join(dataset_t1d_real, by = c("id")) %>%
  
  # remove id variable
  select(-id) %>%
  
  # rename variables
  rename("x" = "recalibration_real", "y" = "recalibration") %>%
  
  # group by real value
  group_by(x) %>%
  
  # create some summary variables
  mutate(
    y_max = quantile(y, probs = c(0.975)),
    y_min = quantile(y, probs = c(0.025))
  ) %>%
  
  # remove grouping
  ungroup() %>%
  
  # create grouping variable
  mutate(group = cut(x, breaks = unique(c(seq(0, 0.1, 0.01), seq(0.1, 0.9, 0.1), seq(0.9, 1, 0.01))))) %>%
  
  # group by groung variable
  group_by(group) %>%
  
  # create summary variables
  mutate(
    x_mean = mean(x),
    ymin = quantile(y, probs = c(0.025)),
    ymin_75 = quantile(y, probs = c(0.125)),
    ymin_50 = quantile(y, probs = c(0.25)),
    ymin_25 = quantile(y, probs = c(0.375)),
    ymax_25 = quantile(y, probs = c(0.625)),
    ymax_50 = quantile(y, probs = c(0.75)),
    ymax_75 = quantile(y, probs = c(0.875)),
    ymax = quantile(y, probs = c(0.975))
  ) %>%
  
  # remove grouping
  ungroup() %>%
  
  # remove variables not needed
  select(-y, -x, -y_max, -y_min, - group) %>%
  
  # keep unique
  distinct()

## custom points for legend
legend_dataset <- cbind(
  x = c(-1, -1, -1, -1),
  y = c(-1, -1, -1, -1),
  alpha = c(0.1, 0.2, 0.3, 0.4)
) %>%
  as.data.frame()


plot_simulation_riley_t1d <- combined_probabilities_t1d %>%
  ggplot() +
  geom_ribbon(aes(x = x_mean, ymin = ymin, ymax = ymax), alpha = 0.1) +
  geom_ribbon(aes(x = x_mean, ymin = ymin_75, ymax = ymax_75), alpha = 0.2) +
  geom_ribbon(aes(x = x_mean, ymin = ymin_50, ymax = ymax_50), alpha = 0.3) +
  geom_ribbon(aes(x = x_mean, ymin = ymin_25, ymax = ymax_25), alpha = 0.4) +
  geom_point(data = legend_dataset, aes(x = x, y = y, alpha = alpha), shape = "square", size = 4) +
  geom_segment(x = 0, xend = 1, y = 0, yend = 1) +
  xlab("Estimated risk from developed model") +
  ylab("Estimated risk from bootstrap models") +
  guides(alpha = guide_legend(title = "Credible interval")) +
  scale_alpha_continuous(labels = c("95%", "75%", "50%", "25%"), range = c(0.1, 0.4)) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
  theme(
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 16),
    title = element_text(size = 17),
    legend.position = "bottom",
    legend.text = element_text(size = 12)) +
  coord_cartesian(ylim = c(0, 1),
                  xlim = c(0, 1))


# Type 2 model
dataset_t2d_real <- data.frame(
  id = 1:nrow(dataset.UNITED_type2_all_genes),
  recalibration_real = predictions_dataset.UNITED_type2_all_genes_new$prob
)

## appending the probabilities together
comparing_t2d_simulation <- NULL

for (i in 1:bootstrap_t2d$iterations) {
  
  comparing_t2d_simulation <- rbind(
    comparing_t2d_simulation,
    bootstrap_t2d$simulations[[i]]$bootstrap_dataset
  )
}

combined_probabilities_t2d <- comparing_t2d_simulation %>%
  
  # join real probability from original model
  left_join(dataset_t2d_real, by = c("id")) %>%
  
  # remove id variable
  select(-id) %>%
  
  # rename variables
  rename("x" = "recalibration_real", "y" = "recalibration") %>%
  
  # group by real value
  group_by(x) %>%
  
  # create some summary variables
  mutate(
    y_max = quantile(y, probs = c(0.975)),
    y_min = quantile(y, probs = c(0.025))
  ) %>%
  
  # remove grouping
  ungroup() %>%
  
  # create grouping variable
  mutate(group = cut(x, breaks = unique(c(seq(0, 0.1, 0.01), seq(0.1, 0.9, 0.1), seq(0.9, 1, 0.01))))) %>%
  
  # group by groung variable
  group_by(group) %>%
  
  # create summary variables
  mutate(
    x_mean = mean(x),
    ymin = quantile(y, probs = c(0.025)),
    ymin_75 = quantile(y, probs = c(0.125)),
    ymin_50 = quantile(y, probs = c(0.25)),
    ymin_25 = quantile(y, probs = c(0.375)),
    ymax_25 = quantile(y, probs = c(0.625)),
    ymax_50 = quantile(y, probs = c(0.75)),
    ymax_75 = quantile(y, probs = c(0.875)),
    ymax = quantile(y, probs = c(0.975))
  ) %>%
  
  # remove grouping
  ungroup() %>%
  
  # remove variables not needed
  select(-y, -x, -y_max, -y_min, - group) %>%
  
  # keep unique
  distinct()

## custom points for legend
legend_dataset <- cbind(
  x = c(-1, -1, -1, -1),
  y = c(-1, -1, -1, -1),
  alpha = c(0.1, 0.2, 0.3, 0.4)
) %>%
  as.data.frame()


plot_simulation_riley_t2d <- combined_probabilities_t2d %>%
  ggplot() +
  geom_ribbon(aes(x = x_mean, ymin = ymin, ymax = ymax), alpha = 0.1) +
  geom_ribbon(aes(x = x_mean, ymin = ymin_75, ymax = ymax_75), alpha = 0.2) +
  geom_ribbon(aes(x = x_mean, ymin = ymin_50, ymax = ymax_50), alpha = 0.3) +
  geom_ribbon(aes(x = x_mean, ymin = ymin_25, ymax = ymax_25), alpha = 0.4) +
  geom_point(data = legend_dataset, aes(x = x, y = y, alpha = alpha), shape = "square", size = 4) +
  geom_segment(x = 0, xend = 1, y = 0, yend = 1) +
  xlab("Estimated risk from developed model") +
  ylab("Estimated risk from bootstrap models") +
  guides(alpha = guide_legend(title = "Credible interval")) +
  scale_alpha_continuous(labels = c("95%", "75%", "50%", "25%"), range = c(0.1, 0.4)) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
  theme(
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 16),
    title = element_text(size = 17),
    legend.position = "bottom",
    legend.text = element_text(size = 12)) +
  coord_cartesian(ylim = c(0, 1),
                  xlim = c(0, 1))


plot_simulation_riley_combined <- wrap_plots(
 
  plot_simulation_riley_t1d +
    ggtitle("Early-insulin-treated"),
  
  plot_simulation_riley_t2d +
    ggtitle("Non-early-insulin-treated"),
  
  nrow = 1
   
) +
  plot_layout(guides = 'collect') &
  theme(
    legend.position = "bottom"
  )


pdf("figures/bootstrap_riley_simulation.pdf", width = 10, height = 5)
plot_simulation_riley_combined
dev.off()
