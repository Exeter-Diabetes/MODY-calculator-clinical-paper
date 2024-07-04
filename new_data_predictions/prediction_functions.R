############################################################################
############################################################################
############################################################################
#
#       This file contains the prediction functions for the MODY models 
#
############################################################################
############################################################################
############################################################################
# 
# # How to make predictions for a new dataset
# 
# # load functions
# source("prediction_functions.R")
# 
# ## load dataset
# ## Format of the dataset - insulin-treated within 6 months
# # #  |  pardm  |  agerec |  hba1c  |  agedx  |    sex     |   bmi   |    C    |    A    |
# # #  |---------|---------|---------|---------|------------|---------|---------|---------|
# # #  | numeric | numeric | numeric | numeric |  numeric   | numeric | numeric | numeric |
# # #  | 1 or 0  |         |         |         | 1 - male   |         | 1 or 0  | 1 or 0  |
# # #  |         |         |         |         | 2 - female |         |         |         |
# # #  |         |         |         |         |            |         |         |         |
# # dataset <- data.frame(cbind(pardm = 1, agerec = 30, hba1c = 7, agedx = 25, sex = 2, bmi = 30, C = 1, A = 0))
# 
# # pardm - parent history of diabetes
# # agerec - age at recruitment
# # hba1c - HbA1c (%)
# # agedx - age at diagnosis
# # sex - sex
# # bmi - BMI
# # insoroha - current treatment on insulin or OHA
# # C - C-peptide positivity or negativity
# # A - Antibody positivity or negativity
# 
# ## Format of the dataset - non-/insulin-treated after 6 months
# # #  |  pardm  |  agerec |  hba1c  |  agedx  |    sex     |   bmi   | insoroha |
# # #  |---------|---------|---------|---------|------------|---------|----------|
# # #  | numeric | numeric | numeric | numeric |  numeric   | numeric | numeric  |
# # #  | 1 or 0  |         |         |         | 1 - male   |         | 1 or 0   |
# # #  |         |         |         |         | 2 - female |         |          |
# # #  |         |         |         |         |            |         |          |
# # dataset <- data.frame(cbind(pardm = 1, agerec = 30, hba1c = 7, agedx = 25, sex = 2, bmi = 30, insoroha = 1))
# 
# 
# ## load data
# # dataset <- ...
# 
# ## load posteriors
# rcs_parms <- readRDS("rcs_parms.rds")
# posterior_samples_T1D <- readRDS("type_1_model_posteriors_single_value.rds")
# 
# # ### create object to use for prediction
# posterior_samples_T1D_obj <- list(post = posterior_samples_T1D$samples)
# class(posterior_samples_T1D_obj) <- "T1D"
# 
# posterior_samples_T2D <- readRDS("type_2_model_posteriors_single_value.rds")
# 
# posterior_samples_T2D_obj <- list(post = posterior_samples_T2D$samples)
# class(posterior_samples_T2D_obj) <- "T2D"
# 
# 
# ## make predictions
# ### For T1D model
# posterior_predictions_T1D <- predict(posterior_samples_T1D_obj, dataset, rcs_parms) %>%
#   apply(., 2, function(x) {
#     data.frame(prob = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
#   }) %>%
#   bind_rows() 
# # 
# # ### For T2D model
# posterior_predictions_T2D <- predict(posterior_samples_T2D_obj, dataset) %>%
#   apply(., 2, function(x) {
#     data.frame(prob = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
#   }) %>%
#   bind_rows() 
# 
# 
#
############################################################################
############################################################################
############################################################################



# load libraries
require(tidyverse)

## Prediction function for MODY T1D

### predict method for 'post' objects
predict.T1D <- function(object, newdata, parms, ...) {
  
  ## check input objects
  stopifnot(class(object) == "T1D")
  stopifnot(is.data.frame(newdata))
  stopifnot(c("pardm", "agerec", "hba1c", "agedx", "sex", "bmi", "C", "A") %in% colnames(newdata))
  # generate variable T
  newdata <- newdata %>%
    mutate(T = ifelse(C == 0 | A == 1, 1, 0))
  x <- select(newdata, one_of(c("pardm", "agerec", "hba1c", "agedx", "sex")))
  xT <- as.matrix(select(newdata, bmi, agedx, pardm, agerec))
  x_spline <- newdata %>%
    select(bmi, agedx, agerec) %>%
    mutate(bmi_spline = as.numeric(rms::rcs(bmi, parms = parms[,"bmi"])[, 2]),
           agedx_spline = as.numeric(rms::rcs(agedx, parms = parms[,"agedx"])[, 2]),
           agerec_spline = as.numeric(rms::rcs(agerec, parms = parms[,"agerec"])[, 2])) %>%
    select(-bmi, -agedx, -agerec) %>%
    set_names(c("bmi_spline", "agedx_spline", "agerec_spline")) %>%
    as.matrix()
  stopifnot(all(map_chr(x, class) == "numeric"))
  stopifnot(is.numeric(newdata$T) | all(is.na(newdata$T)))
  
  ## convert posterior samples to matrix
  post <- as.matrix(do.call(rbind, object$post))
  
  # ## set up data## set up data
  x <- as.matrix(x)
  x <- cbind(int = rep(1, nrow(x)), x)
  betas <- post[, match(c("beta0", paste0("beta[", 1:(ncol(x) - 1), "]")), colnames(post))]
  betas <- x %*% t(betas)
  
  ## do predictions ignoring tests
  preds <- t(t(betas) * post[, "gamma1"])
  preds <- t(t(preds) + post[, "gamma0"])
  preds <- exp(preds) / (1 + exp(preds))
  
  ## all possible combinations of tests
  known_T <- which(!is.na(newdata$T))
  known_nT <- which(is.na(newdata$T))
  
  ## adjust predictions for known T
  if (sum(known_T) > 0) {
    
    # collect values of T and predictions
    T <- newdata$T[!is.na(newdata$T)]
    preds_T <- preds[!is.na(newdata$T),]
    
    # add predictions weighted by variable T
    if (is.null(ncol(preds_T))) {
      pT <- matrix(T, ncol = 1) %*% t(post[, "pMp_Cn_or_Ap"]) + as.vector(1 - T) * preds_T
    } else {
      pT <- matrix(T, ncol = 1) %*% t(post[, "pMp_Cn_or_Ap"]) + matrix(replicate(ifelse(is.null(ncol(preds_T)), 1, ncol(preds_T)), 1 - T), ncol = ifelse(is.null(ncol(preds_T)), 1, ncol(preds_T))) * preds_T
    }
    
    # replace new probability calculations
    preds[!is.na(newdata$T)] <- pT
    
  }
  
  ## adjust predictions for not known T
  if (sum(known_nT) >0) {
    
    # collect values predictions
    preds_nT <- preds[is.na(newdata$T),]
    
    # calculate variable T
    x_all <- cbind(int = rep(1, nrow(x)), xT, x_spline)
    betas <- post[, match(c("beta_t0", paste0("beta_t[", 1:ncol(xT), "]"), paste0("beta_spline[", 1:ncol(x_spline), "]")), colnames(post))]
    betas <- x_all %*% t(betas)
    
    pT_result <- exp(betas) / (1 + exp(betas))
    pT_result <- pT_result[is.na(newdata$T),]
    
    # add predictions weighted by variable T
    if (length(known_nT) == 1) {
      pnT <- pT_result * t(replicate(1, post[, "pMp_Cn_or_Ap"])) + (1 - pT_result)  * preds_nT
    } else {
      pnT <- pT_result * t(replicate(nrow(pT_result), post[, "pMp_Cn_or_Ap"])) + (1 - pT_result)  * preds_nT
    }
    
    # replace new probability calculations
    preds[is.na(newdata$T)] <- pnT
    
  }
  
  ## return posterior predictive samples
  preds <- t(preds)
  colnames(preds) <- NULL
  preds
}


## Prediction function for MODY T2D

### predict method for 'post' objects
predict.T2D <- function(object, newdata, ...) {
  
  #browser()
  
  ## check input objects
  stopifnot(class(object) == "T2D")
  stopifnot(is.data.frame(newdata))
  stopifnot(c("agedx", "bmi", "hba1c", "pardm", "agerec", "insoroha", "sex") %in% colnames(newdata))
  x <- select(newdata, one_of(c("agedx", "bmi", "hba1c", "pardm", "agerec", "insoroha", "sex")))
  stopifnot(all(map_chr(x, class) == "numeric"))
  
  ## convert posterior samples to matrix
  post <- as.matrix(do.call(rbind, object$post))
  
  ## set up data
  x <- as.matrix(x)
  x <- cbind(int = rep(1, nrow(x)), x)
  betas <- post[, match(c("beta0", paste0("beta[", 1:(ncol(x) - 1), "]")), colnames(post))]
  betas <- x %*% t(betas)
  
  ## do predictions ignoring tests
  preds <- t(t(betas) * post[, "gamma1"])
  preds <- t(t(preds) + post[, "gamma0"])
  preds <- exp(preds) / (1 + exp(preds))
  
  ## return posterior predictive samples
  preds <- t(preds)
  colnames(preds) <- NULL
  preds
}
