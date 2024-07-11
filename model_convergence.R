#:--------------------------------------------------------
#
# In this file we plot the model parameter traceplots.
#
#:--------------------------------------------------------

# load libraries
library(tidyverse)
library(ggplot2)
library(patchwork)


# load parameters

type_1_model_posteriors <- readRDS("model_development/type_1_model_posteriors.rds")
type_2_model_posteriors <- readRDS("model_development/type_2_model_posteriors.rds")


# Calculate r hat values
## Turn list of chains into one dataset
type_1_posteriors <- type_1_model_posteriors$samples$chain1 %>%
  rbind(
    type_1_model_posteriors$samples$chain2,
    type_1_model_posteriors$samples$chain3,
    type_1_model_posteriors$samples$chain4 
  ) %>%
  # turn to data.frame
  as.data.frame()

## Calculate r hat
apply(type_1_posteriors, 2, rstan::Rhat)
# R_hat = 1



# Calculate r hat values
## Turn list of chains into one dataset
type_2_posteriors <- type_2_model_posteriors$samples$chain1 %>%
  rbind(
    type_2_model_posteriors$samples$chain2,
    type_2_model_posteriors$samples$chain3,
    type_2_model_posteriors$samples$chain4 
  ) %>%
  # turn to data.frame
  as.data.frame()

## Calculate r hat
apply(type_2_posteriors %>% select(where(is.numeric)), 2, rstan::Rhat)
# R_hat = 1


# labels for parameters
type_1_labels <- list(
  'beta0'="Intercept",
  'beta[1]'="At least one parent affected with diabetes",
  'beta[2]'="Age at recruitment",
  'beta[3]'="HbA1c",
  'beta[4]'="At at diagnosis",
  'beta[5]'="Sex",
  'gamma0'="Gamma 0",
  'gamma1'="Gamma 1",
  'pMp_Cn_or_Ap'="P(M) | C- OR A+"
)
parameter_type_1_labeller <- function(variable,value){
  return(type_1_labels[value])
}


# Trace plots
plot_type_1_model <- type_1_posteriors %>%
  # create variables needed for plotting
  mutate(
    `Chain:` = rep(paste("Chain ", 1:length(type_1_model_posteriors$samples)), each = nrow(type_1_model_posteriors$samples$chain1)),
    `Chain:` = factor(`Chain:`)
  ) %>%
  # gather values into two columns (don't gather chain, used for colour)
  gather(key = "key", value = "value", -`Chain:`) %>%
  # group by chain
  group_by(`Chain:`, key) %>%
  # count the iterations
  mutate(
    Iterations = row_number()
  ) %>%
  # remove grouping
  ungroup() %>%
  # only keep relevant parameters
  filter(key %in% c("beta0", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "pMp_Cn_or_Ap", "gamma0", "gamma1")) %>%
  # order parameters
  mutate(key = factor(key, levels = c("beta0", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "pMp_Cn_or_Ap", "gamma0", "gamma1"))) %>%
  # start plotting
  ggplot(aes(x = Iterations, y = value, colour = `Chain:`)) +
  geom_path(alpha = 0.6) +
  scale_x_continuous(labels = scales::comma) +
  facet_wrap(~key, scales = "free", ncol = 1, labeller=parameter_type_1_labeller) +
  guides(colour = guide_legend(nrow = 1)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title.y = element_blank()
  )
  

pdf("figures/type_1_model_trace_plots.pdf", width = 6, height = 11)
plot_type_1_model
dev.off()


# labels for parameters
type_1_labels <- list(
  'beta0'="Intercept",
  'beta[1]'="At least one parent affected with diabetes",
  'beta[2]'="Age at recruitment",
  'beta[3]'="HbA1c",
  'beta[4]'="At at diagnosis",
  'beta[5]'="Sex",
  'gamma0'="Gamma 0",
  'gamma1'="Gamma 1",
  'pMp_Cn_or_Ap'="P(M) | C- OR A+"
)
# labels for parameters
type_2_labels <- list(
  'beta0'="Intercept",
  'beta[4]'="At least one parent affected with diabetes",
  'beta[5]'="Age at recruitment",
  'beta[3]'="HbA1c",
  'beta[1]'="At at diagnosis",
  'beta[7]'="Sex",
  'beta[6]'="Currently treated with insulin or tablets",
  'beta[2]'="BMI",
  'gamma0'="Gamma 0",
  'gamma1'="Gamma 1"
)
parameter_type_2_labeller <- function(variable,value){
  return(type_2_labels[value])
}

plot_type_2_model <- type_2_posteriors %>%
  # create variables needed for plotting
  mutate(
    `Chain:` = rep(paste("Chain ", 1:length(type_2_model_posteriors$samples)), each = nrow(type_2_model_posteriors$samples$chain1)),
    `Chain:` = factor(`Chain:`)
  ) %>%
  # gather values into two columns (don't gather chain, used for colour)
  gather(key = "key", value = "value", -`Chain:`) %>%
  # group by chain
  group_by(`Chain:`, key) %>%
  # count the iterations
  mutate(
    Iterations = row_number()
  ) %>%
  # remove grouping
  ungroup() %>%
  # order parameters
  mutate(key = factor(key, levels = c("beta0", "beta[4]", "beta[5]", "beta[3]", "beta[1]", "beta[7]", "beta[6]", "beta[2]", "gamma0", "gamma1"))) %>%
  # start plotting
  ggplot(aes(x = Iterations, y = value, colour = `Chain:`)) +
  geom_path(alpha = 0.6) +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(~key, scales = "free", ncol = 1, labeller=parameter_type_2_labeller) +
  guides(colour = guide_legend(nrow = 1)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title.y = element_blank()
  )



pdf("figures/type_2_model_trace_plots.pdf", width = 6, height = 12)
plot_type_2_model
dev.off()

