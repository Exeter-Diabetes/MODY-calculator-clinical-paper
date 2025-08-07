#############################################################################################

#UNITED Discrimination curves

#This contains within it:
# Supplementary Figure 12
# Supplementary Figure 13

#######################################################################################################

# load libraries
library(tidyverse)
library(nimble)
library(pROC)
library(PRROC)


# load functions needed
source("Data/create_data.R")
source("New_Data_Predictions/prediction_functions.R")

# load datasets
## Load population representative dataset
dataset_UNITED_type1 <- create_data(
  dataset = "united t1d", 
  commonmody = FALSE, 
  id = TRUE) %>%
  
  ## if MODY testing missing, change to 0
  mutate(M = ifelse(is.na(M), 0, M))

dataset_UNITED_type2 <- create_data(
  dataset = "united t2d", 
  commonmody = FALSE, 
  id = TRUE)

# load files required
predictions_UNITED_type1_no_T <- readRDS(
  "Model_Predictions/predictions_dataset.UNITED_type1_all_genes_no_T.rds") %>% 
  as.data.frame() %>%
  { rownames(.) <- NULL; . } %>%
  column_to_rownames(var = "id")
predictions_UNITED_type1_no_T <- predictions_UNITED_type1_no_T[as.character(c(dataset_UNITED_type1$id)), ]
predictions_UNITED_type1_with_T <- readRDS(
  "Model_Predictions/predictions_dataset.UNITED_type1_all_genes_with_T.rds") %>% 
  as.data.frame() %>%
  { rownames(.) <- NULL; . } %>%
  column_to_rownames(var = "id")
predictions_UNITED_type1_with_T <- predictions_UNITED_type1_with_T[as.character(c(dataset_UNITED_type1$id)), ]
predictions_UNITED_type2_new <- readRDS(
  "Model_Predictions/predictions_dataset.UNITED_type2_all_genes_new.rds") %>% 
  as.data.frame() %>%
  { rownames(.) <- NULL; . } %>%
  column_to_rownames(var = "id")
predictions_UNITED_type2_new <- predictions_UNITED_type2_new[as.character(c(dataset_UNITED_type2$id)), ]


#:------------------------------------------------------------
# Supplementary Figure 13 -----------------------------------------------------------------
plot_prob_density <- patchwork::wrap_plots(
  
  # T1D models
  patchwork::wrap_plots(
    
    # Clinical features
    patchwork::wrap_plots(
      
      # Density
      dataset_UNITED_type1 %>%
        dplyr::mutate(prob = predictions_UNITED_type1_no_T$prob) %>%
        filter(M == 0) %>%
        select(M, prob) %>%
        rename("Mody" = "M") %>%
        mutate(
          Mody = factor(
            Mody, 
            levels = c(0, 1), 
            labels = c("Non-MODY", "MODY")),
        ) %>%
        ggplot() +
        geom_density(aes(x = prob), fill = "grey") +
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
        ylab("Non-monogenic \n diabetes (n=1,164)"),
      #point
      dataset_UNITED_type1 %>%
        mutate(prob = predictions_UNITED_type1_no_T$prob) %>%
        filter(M == 1) %>%
        select(M, prob) %>%
        rename("Mody" = "M") %>%
        mutate(
          Mody = factor(Mody,
                        levels = c(0, 1), 
                        labels = c("Non-MODY", "MODY")),
        ) %>%
        ggplot() +
        geom_point(aes(x = prob, y=0), 
                   position = position_jitter(height = 0.1, seed = 10)) +
        coord_cartesian(xlim =c(0, 1), ylim = c(-0.15, 0.15)) +
        scale_x_continuous(labels = scales::percent) +
        theme_classic() +
        theme(
          axis.ticks.y = element_blank(), 
          axis.text.y = element_blank(), 
          #axis.title.y = element_blank()
        ) +
        ylab("Monogenic \n diabetes \n (n=10)") +
        xlab("Model probabilities"),
      ncol = 1, nrow = 2, heights = c(4,1)
      
    ),
    
    # Biomarkers
    patchwork::wrap_plots(
      
      # Density
      dataset_UNITED_type1 %>%
        mutate(prob = predictions_UNITED_type1_with_T$prob) %>%
        filter(M == 0) %>%
        select(M, prob) %>%
        rename("Mody" = "M") %>%
        mutate(
          Mody = factor(
            Mody, 
            levels = c(0, 1), 
            labels = c("Non-MODY", "MODY")),
        ) %>%
        ggplot() +
        geom_density(aes(x = prob), fill = "grey") +
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
        ylab("Non-monogenic \n diabetes (n=1,164)"),
      #point
      dataset_UNITED_type1 %>%
        mutate(prob = predictions_UNITED_type1_with_T$prob) %>%
        filter(M == 1) %>%
        select(M, prob) %>%
        rename("Mody" = "M") %>%
        mutate(
          Mody = factor(Mody, 
                        levels = c(0, 1), 
                        labels = c("Non-MODY", "MODY")),
        ) %>%
        ggplot() +
        geom_point(aes(x = prob, y=0), 
                   position = position_jitter(height = 0.1, seed = 24)) +
        coord_cartesian(xlim =c(0, 1), ylim = c(-0.15, 0.15)) +
        scale_x_continuous(labels = scales::percent) +
        theme_classic() +
        theme(
          axis.ticks.y = element_blank(), 
          axis.text.y = element_blank(), 
          #axis.title.y = element_blank()
        ) +
        ylab("Monogenic \n diabetes \n (n=10)") +
        xlab("Model probabilities"),
      ncol = 1, nrow = 2, heights = c(4,1)
      
      
    ),
    
    nrow = 1, ncol = 2
    
  ),
  
  # T2D models
  patchwork::free(patchwork::wrap_plots(
    
    # Density
    dataset_UNITED_type2 %>%
      mutate(prob = predictions_UNITED_type2_new$prob) %>%
      filter(M == 0) %>%
      select(M, prob) %>%
      rename("Mody" = "M") %>%
      mutate(
        Mody = factor(Mody, levels = c(0, 1), labels = c("Non-MODY", "MODY")),
      ) %>%
      ggplot() +
      geom_density(aes(x = prob), fill = "grey") +
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
      ylab("Non-monogenic \n diabetes (n=102)"),
    #point
    dataset_UNITED_type2 %>%
      mutate(prob = predictions_UNITED_type2_new$prob) %>%
      filter(M == 1) %>%
      select(M, prob) %>%
      rename("Mody" = "M") %>%
      mutate(
        Mody = factor(Mody, levels = c(0, 1), labels = c("Non-MODY", "MODY")),
      ) %>%
      ggplot() +
      geom_point(aes(x = prob, y=0), 
                 position = position_jitter(height = 0.1, seed = 14)) +
      coord_cartesian(xlim =c(0, 1), ylim = c(-0.15, 0.15)) +
      scale_x_continuous(labels = scales::percent) +
      theme_classic() +
      theme(
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        #axis.title.y = element_blank()
      ) +
      ylab("Monogenic \n diabetes \n (n=23)") +
      xlab("Model probabilities"),
    ncol = 1, nrow = 2, heights = c(4,1)
    
    
  )
  
  ),
  
  nrow = 2, ncol = 1
  
) + patchwork::plot_annotation(tag_levels = list(c("A.1", "", "A.2", "", "B", ""))) +
  patchwork::plot_layout(
    design = "
    AAAA
    #BB#
    "
  ) &
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 11)
  )


# Making plots
pdf("figures/SF_13_UNITED_prob_density.pdf", width = 13, height = 9)
plot_prob_density
dev.off()

#:------------------------------------------------------------


