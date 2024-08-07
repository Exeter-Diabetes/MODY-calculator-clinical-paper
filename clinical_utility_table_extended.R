##Clinical utility table UNITED -----------------------------------------------------------------
## Extended table --------------------------------------------------------------------------
#load predictions
#load T1D
predictions_dataset.UNITED_type1_new <- readRDS("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper/model_predictions/predictions_dataset.UNITED_type1_with_T.rds")
UNITED_type1 <- cbind(dataset.UNITED_type1_new,predictions_dataset.UNITED_type1_with_T)
#load T2D
predictions_dataset.UNITED_type2_new <- readRDS("~/model_predictions/predictions_dataset.UNITED_type2_new.rds")
UNITED_type2 <- cbind(dataset.UNITED_type2_new,predictions_dataset.UNITED_type2_new)

#Load function to get diagnostic info at different thresholds
calculate_thresholds_diagnostics <- function(response, prediction, unique = FALSE) {
  
  ### threshold number to test
  thresholds <- seq(0, 1, 0.0001)
  
  matrix_thresholds <- matrix(0, nrow = length(thresholds), ncol = 8)
  matrix_thresholds[,1] <- thresholds
  
  preds1 <- data.frame(
    prediction = prediction,
    response = response
  )
  
  for (i in 1:length(thresholds)) {
    # original dataset
    original_data <- preds1
    
    # checking which patients are above the threshold
    interim_above <- preds1 %>%
      mutate(threshold_up = ifelse(prediction > thresholds[i], 1, 0)) %>%
      filter(threshold_up == 1)
    
    interim_below <- preds1 %>%
      mutate(threshold_up = ifelse(prediction > thresholds[i], 1, 0)) %>%
      filter(threshold_up == 0)
    
    # calculating the sensitivity of the MODY cases
    matrix_thresholds[i, 2] <- length(which(interim_above$response == 1))/length(which(original_data$response == 1)) * 100
    
    # calculating number to test
    if (length(which(interim_above$response == 1)) > 0) {
      matrix_thresholds[i, 3] <- ceiling(nrow(interim_above)/length(which(interim_above$response == 1)))
    } else {matrix_thresholds[i, 3] <- as.numeric(NA)}
    
    # Pick-up rate
    matrix_thresholds[i, 4] <- length(which(interim_above$response == 1))/nrow(interim_above) * 100
    
    # calculating specificity
    matrix_thresholds[i, 5] <- length(which(interim_below$response == 0))/length(which(original_data$response == 0)) * 100
    
    # Perc of patients tested
    matrix_thresholds[i, 6] <- (nrow(interim_above)/nrow(original_data))*100
    
    # number of MODY cases missed
    matrix_thresholds[i, 7] <- interim_below %>%
      select(response) %>%
      sum(na.rm = TRUE) %>%
      unlist()
    
    # number of patients tested
    matrix_thresholds[i, 8] <- interim_above %>%
      nrow() %>%
      unlist()
    
  }
  
  ## making the matrix a data.frame()
  matrix_thresholds <- matrix_thresholds %>% 
    as.data.frame() %>%
    set_names(c("Thresholds", "Sensitivity", "NTT", "Pick-up rate", "Specificity", "N-Tested", "Cases Missed", "Patients Tested")) %>%
    arrange(desc(Thresholds))
  
  if (unique == TRUE) {
    ## select unique combinations of sensitivity and NTT (only the first occurance)
    matrix_thresholds <- matrix_thresholds %>%
      slice(which(duplicated(matrix_thresholds %>% select(-Thresholds, -`Pick-up rate`)) == FALSE))
  }
  
  return(matrix_thresholds)
  
}


## Early-insulin-treated -----------------------------------------------------------------
thresholds_UNITED_t1d <- calculate_thresholds_diagnostics(dataset.UNITED_type1_new$M, predictions_dataset.UNITED_type1_with_T$prob)

### 5%
thresholds_UNITED_t1d %>%
  filter(Thresholds == 0.05) %>% arrange(Thresholds) %>% head()


### 10%
thresholds_UNITED_t1d %>%
  filter(Thresholds == 0.1) %>% arrange(Thresholds)  %>% head()



### 20%
thresholds_UNITED_t1d %>%
  filter(Thresholds == 0.2) %>% arrange(Thresholds)  %>% head()



### 25%
thresholds_UNITED_t1d %>%
  filter(Thresholds == 0.25) %>% arrange(Thresholds)  %>% head()



### 30%
thresholds_UNITED_t1d %>%
  filter(Thresholds == 0.3) %>% arrange(Thresholds)  %>% head()



### 36%
thresholds_UNITED_t1d %>%
  filter(Thresholds >= 0.36) %>% arrange(Thresholds)  %>% head()



### 58%
thresholds_UNITED_t1d %>%
  filter(Thresholds >= 0.58) %>% arrange(Thresholds)  %>% head()



### 60%
thresholds_UNITED_t1d %>%
  filter(Thresholds == 0.6) %>% arrange(Thresholds)  %>% head()




## Not-insulin-treated --------------------------------------------------------------------
thresholds_UNITED_t2d <- calculate_thresholds_diagnostics(dataset.UNITED_type2_new$M, predictions_dataset.UNITED_type2_new$prob)

### 5%
thresholds_UNITED_t2d %>%
  filter(Thresholds == 0.05) %>% arrange(Thresholds)  %>% head()


### 10%
thresholds_UNITED_t2d %>%
  filter(Thresholds == 0.1) %>% arrange(Thresholds)  %>% head()



### 20%
thresholds_UNITED_t2d %>%
  filter(Thresholds == 0.2) %>% arrange(Thresholds)  %>% head()



### 25%
thresholds_UNITED_t2d %>%
  filter(Thresholds == 0.25) %>% arrange(Thresholds)  %>% head()


### 30%
thresholds_UNITED_t2d %>%
  filter(Thresholds == 0.3) %>% arrange(Thresholds)  %>% head()



### 36%
thresholds_UNITED_t2d %>%
  filter(Thresholds >= 0.36) %>% arrange(Thresholds)  %>% head()



### 58%
thresholds_UNITED_t2d %>%
  filter(Thresholds >= 0.58) %>% arrange(Thresholds)  %>% head()



### 60%
thresholds_UNITED_t2d %>%
  filter(Thresholds == 0.6) %>% arrange(Thresholds)  %>% head()




