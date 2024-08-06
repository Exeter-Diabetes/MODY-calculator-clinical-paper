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
predictions_dataset.UNITED_type1_with_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_with_T.rds")
predictions_dataset.UNITED_type2_new <- readRDS("model_predictions/predictions_dataset.UNITED_type2_new.rds")


predictions_dataset.UNITED_type1_with_T_full <- readRDS("model_predictions/predictions_dataset.UNITED_type1_with_T_full.rds")
predictions_dataset.UNITED_type2_new_full <- readRDS("model_predictions/predictions_dataset.UNITED_type2_new_full.rds")


# load functions needed
source("data/create_data.R")



## Load population representative dataset
dataset.UNITED_type1 <- create_data(dataset = "united t1d")
# we don't do this to make the plotting easier since these patients don't have the correct biomarker status for testing (cpeptide negative or antibody positive)
# and hence their probabilities are all the same <0.001
  # ## if MODY testing missing, change to 0
  # mutate(M = ifelse(is.na(M), 0, M))    
dataset.UNITED_type2 <- create_data(dataset = "united t2d") %>%
  
  ## if MODY testing missing, change to 0
  mutate(M = ifelse(is.na(M), 0, M))



#:--------------------------------------------

## Calibration plots for UNITED type 1




#
legend_dataset <- cbind(
  x = c(-1, -1, -1, -1),
  y = c(-1, -1, -1, -1),
  alpha = c(0.1, 0.2, 0.3, 0.4)
) %>%
  as.data.frame()

#Define grouping values
#this takes the probabilities in the "probs" column and filters them based on whether their corresponding line in the dataset has a non-missing M
#i.e. only looking at the probabilities for those individuals that have the correct biomarker status (cpeptide positive AND antibody negative)
#that were therefore eligible for testing and will have a prob that isn't (and is greater than) <0.001
grouping_values_UNITED_type1_with_T <- predictions_dataset.UNITED_type1_with_T$prob[which(!is.na(dataset.UNITED_type1$M))] %>% 
  #we then find the 60th and 80th quantile for these individuals
  quantile(probs = c(0.2, 0.4, 0.6, 0.8))

predictions_dataset.UNITED_type1_with_T$prob[which(!is.na(dataset.UNITED_type1$M))] %>% 
  #we then find the 60th and 80th quantile for these individuals
  quantile(probs = c(0.2, 0.4, 0.6, 0.8))

predictions_dataset.UNITED_type1_with_T$prob[which(!is.na(dataset.UNITED_type1$M))] %>% 
  #we then find the 60th and 80th quantile for these individuals
  quantile(probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))

plot_calibration_UNITED_type1_with_T <- predictions_dataset.UNITED_type1_with_T_full[,which(!is.na(dataset.UNITED_type1$M))] %>%
  as.data.frame() %>%
  #transpose data - swop rows and columns
  t() %>%
  as.data.frame() %>%
  #gather() does the same as pivot_longer
  #i.e. all column names become categories under the new "key" column 
  #and all the values that went in the respective columns get moved to the corresponding space in the "values" column
  gather() %>%
  #bind two columns in addition to the "key" and "value" columns
  cbind(
    #group denotes the predicted probabilities for those with MODY testing (& correct biomarker status)
    group = predictions_dataset.UNITED_type1_with_T$prob[which(!is.na(dataset.UNITED_type1$M))],
    #M denotes the MODY status (1= has MODY gene; 0 = doesn't have MODY gene)
    M = dataset.UNITED_type1$M[!is.na(dataset.UNITED_type1$M)]
  ) %>%
  #rewrite the group variable to groups (1-5) 
  #corresponding to the predicted probabilities orientation around the quantiles defined in grouping_values
  mutate(
    #e.g.if the predicted probs are less than the 1st element in grouping_values (quantile), then change its value to 1
    group = ifelse(group < grouping_values_UNITED_type1_with_T[1], 1,
                   ifelse(group < grouping_values_UNITED_type1_with_T[2], 2,
                          ifelse(group < grouping_values_UNITED_type1_with_T[3], 3,
                                 ifelse(group < grouping_values_UNITED_type1_with_T[4], 4, 5))))
  ) %>%
  #now we want to group by these 5 defined groups (or defined predicted probability ranges)
  group_by(group) %>%
  #now we want to define new variables for each group
  mutate(
    #x which is the 50th quantile of the probs in value for each group (median probability)
    x = quantile(value, probs = c(0.5)),
    #y which is the proportion of MODY cases in each group
    y = sum(M)/n(),
    #xmin which is the 2.5th quantile contributing to the lower limit of the 95% CI bounds
    #xmin and max provide 95% CI
    xmin = quantile(value, probs = c(0.025)),
    #xmin75 provides the 12.5th quantile contributing to the lower limit of the 75% CI bounds
    xmin_75 = quantile(value, probs = c(0.125)),
    #xmin50 which is the 25th quantile contributing to the lower limit of the 50% CI bounds
    xmin_50 = quantile(value, probs = c(0.25)),
    #xmin25 provides the 37.5th quantile contributing to the lower limit of the 25% CI bounds
    xmin_25 = quantile(value, probs = c(0.375)),
    #xmax25 provides the 62.5th quantile contributing to the upper limit of the 25% CI bounds
    xmax_25 = quantile(value, probs = c(0.625)),
    #xmax50 which is the 75th quantile contributing to the upper limit of the 50% CI bounds
    xmax_50 = quantile(value, probs = c(0.75)),
    #xmax75 which is the 87.5th quantile contributing to the upper limit of the 75% CI bounds
    xmax_75 = quantile(value, probs = c(0.875)),
    #xmax which is the 97.5th quantile contributing to the upper limit of the 95% CI bounds
    xmax = quantile(value, probs = c(0.975))
  ) %>%
  #ungroup by group
  ungroup() %>%
  #we now have specific x,y, xmin etc for each group (quintile)
  #i.e. everytime see a group 4 will have same x, y, xmin etc
  #as these are the key values needed to plot, we don't need the rest and therefore deselect/drop them
  select(-key, -value, -group, -M) %>%
  #currently have 80000000 rows of repeating rows
  #use distinct to keep only unique rows
  distinct() %>%
  #now we can start plotting
  #x is the median predicted probability
  #y is the proportion of MODY in that quintile
  ggplot(aes(x = x, y = y)) +
  #PLOTTING
  #plot y=x reference line
  geom_abline(aes(intercept = 0, slope = 1)) +
  #this defines the 95% CI
  geom_ribbon(aes(xmin = xmin, xmax = xmax), alpha = 0.1) +
  #this defines the 75% CI
  geom_ribbon(aes(xmin = xmin_75, xmax = xmax_75), alpha = 0.2) +
  #the 50% CI
  geom_ribbon(aes(xmin = xmin_50, xmax = xmax_50), alpha = 0.3) +
  #the 25% CI
  geom_ribbon(aes(xmin = xmin_25, xmax = xmax_25), alpha = 0.4) +
  #plot a linear line of best fit
  #formula = y ~x is the default smoothing function
  geom_smooth(method = 'lm', formula = y ~ x,
              aes(fill = after_scale(color)), alpha = 0) +
  #plot x and y points as triangles
  geom_point(shape = "triangle") +
  #AXES FORMATTING
  #this ensures the ribbons are plotted correctly
  #limit x axis between 0 and 1
  xlim(0, 1) +
  #limit y axis between 0 and 1
  ylim(0, 1) +
  #set the coordinated of the cartesian plane (axes)
  #this allows us to zoom into the area of interest without messing with the ribbon plotting
  coord_cartesian(ylim =c(0, 0.3), xlim =c(0, 1)) +
  #rewrite the x and y axes scales to percent for the purpose of labelling
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  #Provide axis labels
  xlab("Model predictions") +
  ylab("Observed probability") +
  #LEGEND FORMATTING
  #this provides the basis for the legend mapping alpha values to credible intervals
  #points with assigned alphas are plotted outside of the axis limits
  #this is to have a reference for ggplot to form the legend, but without it appearing on the plot
  geom_point(data = legend_dataset, aes(x = x, y = y, alpha = alpha), shape = "square", size = 4) +
  #provides a legend for the scaled alpha values
  #this says the title of the legend guide for the alpha values, which are mapped to the credible intervals
  #guides(alpha = guide_legend(title = "Credible interval")) +
  #this provides the name of each element in the alpha-credible interval scale that we want in the legend
  #and the corresponding alpha ranges that they fall within
  #scale_alpha_continuous(labels = c("95%", "75%", "50%", "25%"), range = c(0.1, 0.4)) +
  scale_alpha_continuous(name = "Credible interval", labels = c("95%", "75%", "50%", "25%"), range = c(0.1, 0.4)) +
  #THEME FORMATTING
  theme_light() +
  theme(
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 16),
    title = element_text(size = 17),
    legend.position = "bottom",
    legend.text = element_text(size = 12))

plot(plot_calibration_UNITED_type1_with_T)

#Define grouping values
#this takes the probabilities in the "probs" column and filters them based on whether their corresponding line in the dataset has a non-missing M
#i.e. only looking at the probabilities for those individuals weren't tested for MODY
#in this case this is irrelevant as there are none and we altered any missing to 0 
#therefore looking for quantiles on all predicted probs
grouping_values_UNITED_type2_new <- predictions_dataset.UNITED_type2_new$prob[which(!is.na(dataset.UNITED_type2$M))] %>% 
  #we then find the 20th and 80th quantile for these individuals
  quantile(probs = c(0.2, 0.8))



plot_calibration_UNITED_type2_new <- predictions_dataset.UNITED_type2_new_full %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  gather() %>%
  cbind(
    group = predictions_dataset.UNITED_type2_new$prob,
    M = dataset.UNITED_type2$M
  ) %>%
  mutate(
    group = ifelse(group < grouping_values_UNITED_type2_new[1], 1,
                   ifelse(group < grouping_values_UNITED_type2_new[2], 2,
                          ifelse(group < grouping_values_UNITED_type2_new[3], 3,
                                 ifelse(group < grouping_values_UNITED_type2_new[4], 4, 5))))
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
  coord_cartesian(ylim =c(0, 0.6), xlim =c(0, 1)) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  theme(
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 16),
    title = element_text(size = 17),
    legend.position = "bottom",
    legend.text = element_text(size = 12))





plot_calibration <- patchwork::wrap_plots(
  
  plot_calibration_UNITED_type1_with_T,
  
  plot_calibration_UNITED_type2_new,
  
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
