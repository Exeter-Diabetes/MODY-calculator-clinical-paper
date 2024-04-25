#:--------------------------------------------------------
#   
# In this file we check the probabilities for CPRD individuals
#
#:--------------------------------------------------------

# load libraries
library(tidyverse)
library(ggplot2)
library(patchwork)


# load data
probabilities_under_30 <- read.delim("data/mody_probabilities_under_30s.txt")
probabilities_under_35 <- read.delim("data/mody_probabilities_under_35s.txt")


## formatted probabilities
probabilities <- probabilities_under_30 %>%
  
  as.data.frame() %>%
  
  ## remove the column which the model used
  select(-which_equation) %>%
  
  ## rename variables
  rename("Old Calculator" = "mody_adj_prob_fh0", "New Calculator" = "pedro_prob") %>%
  
  ## gather in long format
  gather("key", "value") %>%
  
  ## turn key into factor
  mutate(key = factor(key, levels = c("Old Calculator", "New Calculator"))) %>%
  
  ## divide values by 100 to decimals
  mutate(value = value/100)




plot_cprd_probabilities <- patchwork::wrap_plots(
  
  probabilities %>%
    ggplot() +
    geom_vline(xintercept = 0.036, colour = "black") +
    geom_boxplot(aes(x = value, y = key, colour = key)) +
    geom_point(
      data = probabilities %>%
        group_by(key) %>%
        mutate(mean = mean(value)) %>%
        ungroup() %>%
        select(-value) %>%
        unique(),
      aes(x = mean, y = key), shape = 18, size = 4
    ) +
    scale_colour_manual(values = c("#1B9E77", "#D95F02")) +
    coord_cartesian(xlim = c(0, 1)) +
    scale_x_continuous(labels = scales::percent, breaks = c(0, 0.01, 0.025, 0.05, 0.1, 0.2, 0.4, 1)) +
    coord_trans(x = "log2") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 22),
      axis.title = element_blank(),
      axis.text = element_text(size = 22),
      legend.position = "none",
      axis.line = element_line(colour = "black"),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks.y = element_blank()
    ),
  
  probabilities %>%
    ggplot() +
    geom_vline(xintercept = 0.036, colour = "black") +
    geom_boxplot(aes(x = value, y = key, colour = key)) +
    geom_point(
      data = probabilities %>%
        group_by(key) %>%
        mutate(mean = mean(value)) %>%
        ungroup() %>%
        select(-value) %>%
        unique(),
      aes(x = mean, y = key), shape = 18, size = 4
    ) +
    scale_colour_manual(values = c("#1B9E77", "#D95F02")) +
    coord_cartesian(xlim = c(0, 1)) +
    scale_x_continuous(labels = scales::percent) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 22),
      axis.title = element_blank(),
      axis.text = element_text(size = 22),
      legend.position = "none",
      axis.line = element_line(colour = "black"),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks.y = element_blank()
    ),
  ncol = 1
)


