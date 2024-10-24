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
  rename("Old Calculator" = "old_prob", "New Calculator" = "new_prob") %>%
  
  ## gather in long format
  gather("key", "value") %>%
  
  ## turn key into factor
  mutate(key = factor(key, levels = c("Old Calculator", "New Calculator")))
  
  ## divide values by 100 to decimals
  # mutate(value = value/100)


plot_cprd_probabilities_density <- patchwork::wrap_plots(
  
  # the first part of the plot
  probabilities %>%
    filter(key == "New Calculator") %>%
    ggplot() +
    geom_density(aes(x = value), fill = "grey") +
    geom_vline(
      data = probabilities %>%
        filter(key == "New Calculator") %>%
        group_by(key) %>%
        mutate(mean = mean(value)) %>%
        ungroup() %>%
        select(-value) %>%
        unique(),
      aes(xintercept = mean), linetype = "dashed"
    ) +
    geom_vline(xintercept = 0.036, colour = "black") +
    coord_cartesian(xlim = c(0, 0.115), expand = TRUE) +
    scale_x_continuous(labels = scales::percent, breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 22),
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none",
      axis.line = element_line(colour = "black"),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = unit(c(0,0,0,0), "cm")
    ),
  
  # spacing
  plot_spacer(),
  
  # the second part of the plot
  probabilities %>%
    filter(key == "New Calculator") %>%
    ggplot() +
    geom_density(aes(x = value), fill = "grey") +
    geom_vline(
      data = probabilities %>%
        filter(key == "New Calculator") %>%
        group_by(key) %>%
        mutate(mean = mean(value)) %>%
        ungroup() %>%
        select(-value) %>%
        unique(),
      aes(xintercept = mean), linetype = "dashed"
    ) +
    geom_vline(xintercept = 0.036, colour = "black") +
    coord_cartesian(xlim = c(0.16, 1), expand = TRUE) +
    scale_x_continuous(labels = scales::percent, breaks = c(0.12, 0.25, 0.5, 0.75, 1)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 22),
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none",
      axis.line = element_line(colour = "black"),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      plot.margin = unit(c(0,0.2,0,0), "cm")
    )
  
) + plot_layout(widths = c(4, -0.17, 4))



plot_cprd_probabilities_boxplot <- patchwork::wrap_plots(
  
  probabilities %>%
    ggplot() +
    geom_boxplot(aes(x = value, y = key, colour = key)) +
    geom_point(
      data = probabilities %>%
        group_by(key) %>%
        mutate(mean = mean(value)) %>%
        ungroup() %>%
        select(-value) %>%
        unique(),
      aes(x = mean, y = key, colour = key), shape = "circle", size = 5
    ) +
    geom_vline(xintercept = 0.036, colour = "black") +
    scale_colour_manual(values = c("#1B9E77", "#D95F02")) +
    coord_cartesian(xlim = c(0, 0.115), expand = TRUE) +
    scale_x_continuous(labels = scales::percent, breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 22),
      axis.title = element_blank(),
      axis.text.y = element_text(size = 22),
      legend.position = "none",
      axis.line = element_line(colour = "black"),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = unit(c(0,1.5,0,0), "cm")
    ),
  
  probabilities %>%
    ggplot() +
    geom_boxplot(aes(x = value, y = key, colour = key)) +
    geom_point(
      data = probabilities %>%
        group_by(key) %>%
        mutate(mean = mean(value)) %>%
        ungroup() %>%
        select(-value) %>%
        unique(),
      aes(x = mean, y = key, colour = key), shape = "circle", size = 5
    ) +
    geom_vline(xintercept = 0.036, colour = "black") +
    scale_colour_manual(values = c("#1B9E77", "#D95F02")) +
    coord_cartesian(xlim = c(0.16, 1), expand = TRUE) +
    scale_x_continuous(labels = scales::percent, breaks = c(0.125, 0.25, 0.5, 0.75, 1)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 22),
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none",
      axis.line = element_line(colour = "black"),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      plot.margin = unit(c(0,0.5,0,0), "cm")
    ),
  
  nrow = 1
)



#:-------------------------------------------------------------
# Making plots
pdf("figures/cprd_analysis_density.pdf", width = 11, height = 3)
plot_cprd_probabilities_density
dev.off()



