#:--------------------------------------------------------
#   
# In this file we check convergence in model parameters for the calculator.
#
#:--------------------------------------------------------

# load libraries
library(tidyverse)
library(nimble)
library(patchwork)

# load files required
type_1_model_posteriors <- readRDS("model_development/type_1_model_posteriors.rds")
type_2_model_posteriors <- readRDS("model_development/type_2_model_posteriors.rds")
type_1_old_model_posteriors <- readRDS("model_development/type_1_old_model_posteriors.rds")
type_2_old_model_posteriors <- readRDS("model_development/type_2_old_model_posteriors.rds")


# load function needed
long_list <- function(data) {
  
  # check number of chains in data
  num_chain <- length(data)
  
  # this file will contain the output
  output <- NULL
  
  # iterate through each chain
  for (i in 1:num_chain) {
    
    # append below the output file
    output <- rbind(
      output,
      
      # original data
      data[[i]] %>%
        
        # turn to matrix
        as.matrix() %>%
        
        # turn to data.frame
        as.data.frame() %>%
        
        # turn everything into two columns
        gather() %>%
        
        # add columns needed
        mutate(
          
          # Chain number
          Chain = paste0("Chain ", i),
          
          # number of iteration
          iteration = rep(1:nrow(data[[i]] %>% as.matrix()), ncol(data[[i]] %>% as.matrix()))
          
        )
    )
    
  }
  
  # return the output file
  return(output %>% as.data.frame())
  
}

#:---------------------------------------------------------------------------------------------------
#:---------------------------------------------------------------------------------------------------

# trace plot for old model type 1
plot_old_type_1_values <- long_list(type_1_old_model_posteriors$samples) %>%
  
  # modify variable to include the right labels
  mutate(key = factor(
    
    ## column name
    key, 
    
    ## select the levels for the variable
    levels = c("beta0", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]"),
    
    # select the labels for the variable levels
    labels = c("Intercept", "At least one parent affected with diabetes", "Age at recruitment (years)", "Hba1c (%)", "Age at diagnosis (years)", "Sex (baseline Male)"))
  ) %>%
  
  # start ggplot()
  ggplot() +
  
  # plot the path trace
  geom_path(aes(x = iteration, y = value, colour = Chain)) +
  
  # divide the plot into rows of each parameter
  facet_wrap(~key, ncol = 1, scales = "free") +
  
  # change the x axis to have commas
  scale_x_continuous("Iterations", label=scales::comma) +
  
  # choose the theme for the plot
  theme_light() +
  
  # choose the palette
  scale_color_brewer(palette="Set1") +
  
  # specific changes to plot details
  theme(
    
    ## no legend
    legend.position = "none",
    
    ## no title for y axis
    axis.title.y = element_blank(),
    
    ## title size for x axis
    axis.title.x = element_text(size = 10),
    
    ## text size for both axis
    axis.text = element_text(size = 10),
    
    ## text size for the title strips for facet_wrap 
    strip.text.x = element_text(size = 15)
  )


# trace plot for model type 2
plot_old_type_2_values <- long_list(type_2_old_model_posteriors$samples) %>%
  
  # modify variable to include the right labels
  mutate(key = factor(
    
    ## column name
    key,
    
    ## select the levels for the variable
    levels = c("beta0", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "beta[6]", "beta[7]"),
    
    # select the labels for the variable levels
    labels = c("Intercept", "Age at diagnosis (years)", "BMI (kg/m2)", "Hba1c (%)", "At least one parent affected with diabetes", 
               "Age at recruitment (years)", "Insulin or OHA tablet", "Sex (baseline Male)"))
  ) %>%
  
  # start ggplot()
  ggplot() +
  
  # plot the path trace
  geom_path(aes(x = iteration, y = value, colour = Chain)) +
  
  # divide the plot into rows of each parameter
  facet_wrap(~key, ncol = 1, scales = "free") +
  
  # change the x axis to have commas
  scale_x_continuous("Iterations", label=scales::comma) +
  
  # choose the theme for the plot
  theme_light() +
  
  # choose the palette
  scale_color_brewer(palette="Set1") +
  
  # specific changes to plot details
  theme(
    
    ## no legend
    legend.position = "none",
    
    ## no title for y axis
    axis.title.y = element_blank(),
    
    ## title size for x axis
    axis.title.x = element_text(size = 10),
    
    ## text size for both axis
    axis.text = element_text(size = 10),
    
    ## text size for the title strips for facet_wrap 
    strip.text.x = element_text(size = 15)
  )

# combine both plots into one
plot_combined_old_calculator <- patchwork::wrap_plots(
  
  # plot for type 1 model
  plot_old_type_1_values,
  
  # plot for type 2 model
  plot_old_type_2_values,
  
  # new plot layout
  ncol = 2, nrow = 1
) + 
  
  # add A/B labels
  patchwork::plot_annotation(tag_levels = 'A')



#:---------------------------------------------------------------------------------------------------
#:---------------------------------------------------------------------------------------------------


# trace plot for model type 1
plot_type_1_values <- long_list(type_1_model_posteriors$samples) %>%
  
  # modify variable to include the right labels
  mutate(key = factor(
    
    ## column name
    key, 
    
    ## select the levels for the variable
    levels = c("beta0", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "gamma0" , "gamma1", "pMp_Cn_or_Ap",
               "beta_t0", "beta_t[1]", "beta_spline[1]", "beta_t[2]", "beta_spline[2]", "beta_t[4]", "beta_spline[3]", "beta_t[3]"),
    
    # select the labels for the variable levels
    labels = c("Intercept", "At least one parent affected with diabetes", "Age at recruitment (years)", "Hba1c (%)", "Age at diagnosis (years)", "Sex (baseline Male)", "Gamma 0", "Gamma 1", "P(MODY+|C- OR A+)",
               "T - Intercept", "T - BMI (kg/m2)", "T - Spline(BMI)", "T - Age at diagnosis (years)", "T - Spline(Age at diagnosis)", "T - Age at recruitment (years)", "T - Spline(Age at recruitment)", "T - At least one parent affected with diabetes"))
  ) %>%
  
  # start ggplot()
  ggplot() +
  
  # plot the path trace
  geom_path(aes(x = iteration, y = value, colour = Chain)) +
  
  # divide the plot into rows of each parameter
  facet_wrap(~key, ncol = 1, scales = "free") +
  
  # change the x axis to have commas
  scale_x_continuous("Iterations", label=scales::comma) +
  
  # choose the theme for the plot
  theme_light() +
  
  # choose the palette
  scale_color_brewer(palette="Set1") +
  
  # specific changes to plot details
  theme(
    
    ## no legend
    legend.position = "none",
    
    ## no title for y axis
    axis.title.y = element_blank(),
    
    ## title size for x axis
    axis.title.x = element_text(size = 10),
    
    ## text size for both axis
    axis.text = element_text(size = 10),
    
    ## text size for the title strips for facet_wrap 
    strip.text.x = element_text(size = 15)
  )


# trace plot for model type 2
plot_type_2_values <- long_list(type_2_model_posteriors$samples) %>%
  
  # modify variable to include the right labels
  mutate(key = factor(
    
    ## column name
    key,
    
    ## select the levels for the variable
    levels = c("beta0", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "beta[6]", "beta[7]", "gamma0" , "gamma1"),
    
    # select the labels for the variable levels
    labels = c("Intercept", "Age at diagnosis (years)", "BMI (kg/m2)", "Hba1c (%)", "At least one parent affected with diabetes", 
               "Age at recruitment (years)", "Insulin or OHA tablet", "Sex (baseline Male)", "Gamma 0", "Gamma 1"))
  ) %>%
  
  # start ggplot()
  ggplot() +
  
  # plot the path trace
  geom_path(aes(x = iteration, y = value, colour = Chain)) +
  
  # divide the plot into rows of each parameter
  facet_wrap(~key, ncol = 1, scales = "free") +
  
  # change the x axis to have commas
  scale_x_continuous("Iterations", label=scales::comma) +
  
  # choose the theme for the plot
  theme_light() +
  
  # choose the palette
  scale_color_brewer(palette="Set1") +
  
  # specific changes to plot details
  theme(
    
    ## no legend
    legend.position = "none",
    
    ## no title for y axis
    axis.title.y = element_blank(),
    
    ## title size for x axis
    axis.title.x = element_text(size = 10),
    
    ## text size for both axis
    axis.text = element_text(size = 10),
    
    ## text size for the title strips for facet_wrap 
    strip.text.x = element_text(size = 15)
  )

# combine both plots into one
plot_combined_new_calculator <- patchwork::wrap_plots(
  
  # plot for type 1 model
  plot_type_1_values,
  
  # plot for type 2 model
  plot_type_2_values,
  
  # new plot layout
  ncol = 2, nrow = 1
) + 
  
  # add A/B labels
  patchwork::plot_annotation(tag_levels = 'A')


