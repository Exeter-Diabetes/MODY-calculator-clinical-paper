# MYDIABETES NEW MODELS WHITE ETHNICITY -----------------------------------------------
#SET WORKING DIRECTORY -----------------------------------------------------------------------------------
setwd("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper/new_data_predictions")
#load libraries ------------------------------------------------------------------------------------
library(tidyverse)
library(writexl)
library(pROC)

#readin datasets --------------------------------------------------------------------------------------
#load T1D
load("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper/new_data_predictions/MY_T1D.RData")
#load T2D 
load("~/PhD/CLINICAL MODY/Code/MODY-calculator-clinical-paper/new_data_predictions/MY_T2D.RData")

# LIMIT TO ONLY WHITE ETHNICITY
MY_T1D <- MY_T1D %>%
  filter(ethnicity == "White")

MY_T2D <- MY_T2D %>%
  filter(ethnicity == "White")

#DATASETS PREP --------------------------------------------------------------------------------------
#T1D data
MY_T1D <- MY_T1D %>%
  rename(
    hba1c_mmol = hba1c
  )

MY_T1D <- MY_T1D %>%
  rename(
    hba1c = hba1c_perc, 
    bmi = BMI
  )
#T2D data
MY_T2D <- MY_T2D %>%
  rename(
    hba1c_mmol = hba1c
  )

MY_T2D <- MY_T2D %>%
  rename(
    hba1c = hba1c_perc, 
    bmi = BMI
  )

#PREDICTIONS PREP --------------------------------------------------------------------------------
# load functions
source("prediction_functions.R")

# ## load posteriors
rcs_parms <- readRDS("rcs_parms.rds")

##T1D/EARLY-INSULIN-TREATED MODEL -------------------------------------------------------------------
posterior_samples_T1D <- readRDS("type_1_model_posteriors_thin_100.rds")
#
# ### create object to use for prediction
posterior_samples_T1D_obj <- list(post = posterior_samples_T1D$samples)
class(posterior_samples_T1D_obj) <- "T1D"

## make predictions
posterior_predictions_T1D <- predict(posterior_samples_T1D_obj, MY_T1D, rcs_parms) %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, prob = 0.025), UCI = quantile(x, prob = 0.975))
  }) %>%
  bind_rows()

#MERGE TO MY_T1D
MY_T1D <- cbind(MY_T1D, posterior_predictions_T1D)

##T2D/NOT-EARLY-INSULIN-TREATED MODEL -----------------------------------------------------------------
posterior_samples_T2D <- readRDS("type_2_model_posteriors_thin_100.rds")
# 
posterior_samples_T2D_obj <- list(post = posterior_samples_T2D$samples)
class(posterior_samples_T2D_obj) <- "T2D"

## make predictions
posterior_predictions_T2D <- predict(posterior_samples_T2D_obj, MY_T2D) %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, prob = 0.025), UCI = quantile(x, prob = 0.975))
  }) %>%
  bind_rows()

#MERGE TO MY_T2D
MY_T2D <- cbind(MY_T2D, posterior_predictions_T2D)

## Early-insulin-treated model analysis -----------------------------------------------------------------------
MY_T1D %>%
  group_by(M) %>%
  summarise(mean = mean(prob),
            meanLCI = mean(LCI),
            meanUCI = mean(UCI),
            n= n())


##BOXPLOTS ----------------------------------------------------------------------------------------
ggplot(MY_T1D, aes(x = factor(M, levels = c(0, 1), labels = c("Negative", "Positive")), y = prob)) +
  geom_boxplot() +
  labs(x = "MODY Status", y = "Probability") +
  theme_bw()

##ROC -----------------------------------------------------------------------------------------

roc_t1d_model <- roc(MY_T1D$M, MY_T1D$prob, plot = TRUE, print.thres = "best", print.auc = TRUE)
model_t1d_pr <- coords(roc_model, x = "best", ret=c("threshold", "specificity", "sensitivity", "accuracy", "precision", "recall", "ppv", "npv"), transpose = FALSE)
#model_pr <- model_pr[c(-1,-3),]


sum_t1dtable$threshold <- model_pr$threshold
sum_t1dtable$ROCAUC <- roc_model$auc
sum_t1dtable$accuracy <- model_pr$accuracy
sum_t1dtable$sensitivity <- model_pr$sensitivity
sum_t1dtable$specificity <- model_pr$specificity
sum_t1dtable$PPV <- model_pr$ppv
sum_t1dtable$NPV <- model_pr$npv

print(roc_t1d_model)

MY_T1D_mydiabetes_ROC_table <- as.data.frame(sum_t1dtable)
write_xlsx(MY_T1D_mydiabetes_ROC_table,"MY_T1D_mydiabetes_ROC_table.xlsx")

## Calibration ---------------------------------------------------------------------------------------
brks_t1d <- unique(quantile(MY_T1D_WHITE$prob, prob = seq(0, 1, by = 0.2), na.rm = TRUE))
brks_t1d[1] <- 0.99 * brks_t1d[1]
brks_t1d[length(brks_t1d)] <- 1.1 * brks_t1d[length(brks_t1d)]
dec_t1d <- cut(
  MY_T1D_WHITE$prob, 
  breaks = brks_t1d,
  include_lowest = TRUE
)
  
dec_t1d <- data.frame(y = MY_T1D_WHITE$M, pred = MY_T1D_WHITE$prob, dec = dec_t1d) %>%
  group_by(dec) %>%
  mutate(prob_obs = sum(y) / n(), 
          obs = sum(y),
          n_group = n(),
          mnpred = mean(pred),
          lower = lapply(sum(y), prop.test, n = n()), 
          upper = sapply(lower, function(x) x$conf.int[2]), 
          lower = sapply(lower, function(x) x$conf.int[1]))
  
  ## plot 
ggplot(dec_t1d, aes(x = mnpred, y = prob_obs)) +
  geom_point() +
  xlab("Mean predicted probability in each decile") +
  ylab("Proportion of MODY in each decile") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  #ylim(c(0, 1)) + xlim(c(0, 1)) +
  geom_errorbar(aes(ymin = lower, ymax = upper))

##Thresholds -----------------------------------------------------------------------------
### 5%
MY_T1D_WHITE_5PERC <- MY_T1D %>%
  summarise(totalover = sum(prob >= 0.05),
            ncasespickedup = sum(prob >= 0.05 & M ==1),
            PPV = (sum(prob >= 0.05 & M ==1)/sum(prob >= 0.05))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.05 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.05 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.05 & M == 0)/sum(prob < 0.05))*100,
            Sensitivity = (sum(prob >= 0.05 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.05 & M == 0)/sum(M==0))*100)
### 10%
MY_T1D_WHITE_10PERC <- MY_T1D %>%
  summarise(totalover = sum(prob >= 0.1),
            ncasespickedup = sum(prob >= 0.1 & M ==1),
            PPV = (sum(prob >= 0.1 & M ==1)/sum(prob >= 0.1))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.1 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.1 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.1 & M == 0)/sum(prob < 0.1))*100,
            Sensitivity = (sum(prob >= 0.1 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.1 & M == 0)/sum(M==0))*100)

### 20%
MY_T1D_WHITE_20PERC <- MY_T1D %>%
  summarise(totalover = sum(prob >= 0.2),
            ncasespickedup = sum(prob >= 0.2 & M ==1),
            PPV = (sum(prob >= 0.2 & M ==1)/sum(prob >= 0.2))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.2 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.2 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.2 & M == 0)/sum(prob < 0.2))*100,
            Sensitivity = (sum(prob >= 0.2 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.2 & M == 0)/sum(M==0))*100)

MY_T1D_WHITE_THRESHOLDS <- rbind(MY_T1D_WHITE_5PERC, MY_T1D_WHITE_10PERC, MY_T1D_WHITE_20PERC)
write_xlsx(MY_T1D_WHITE_THRESHOLDS,"MY_T1D_WHITE_THRESHOLDS_table.xlsx")

## Early-insulin-treated model analysis -----------------------------------------------------------------------
MY_T2D %>%
  group_by(M) %>%
  summarise(mean = mean(prob),
            meanLCI = mean(LCI),
            meanUCI = mean(UCI),
            n= n())


##BOXPLOTS ----------------------------------------------------------------------------------------
ggplot(MY_T2D, aes(x = factor(M, levels = c(0, 1), labels = c("Negative", "Positive")), y = prob)) +
  geom_boxplot() +
  labs(x = "MODY Status", y = "Probability") +
  theme_bw()

##ROC -----------------------------------------------------------------------------------------

roc_T2D_model <- roc(MY_T2D$M, MY_T2D$prob, plot = TRUE, print.thres = "best", print.auc = TRUE)
model_T2D_pr <- coords(roc_model, x = "best", ret=c("threshold", "specificity", "sensitivity", "accuracy", "precision", "recall", "ppv", "npv"), transpose = FALSE)
#model_pr <- model_pr[c(-1,-3),]


sum_T2Dtable$threshold <- model_pr$threshold
sum_T2Dtable$ROCAUC <- roc_model$auc
sum_T2Dtable$accuracy <- model_pr$accuracy
sum_T2Dtable$sensitivity <- model_pr$sensitivity
sum_T2Dtable$specificity <- model_pr$specificity
sum_T2Dtable$PPV <- model_pr$ppv
sum_T2Dtable$NPV <- model_pr$npv

print(roc_T2D_model)

MY_T2D_mydiabetes_ROC_table <- as.data.frame(sum_T2Dtable)
write_xlsx(MY_T2D_mydiabetes_ROC_table,"MY_T2D_mydiabetes_ROC_table.xlsx")

## Calibration ---------------------------------------------------------------------------------------
brks_T2D <- unique(quantile(MY_T2D_WHITE$prob, prob = seq(0, 1, by = 0.2), na.rm = TRUE))
brks_T2D[1] <- 0.99 * brks_T2D[1]
brks_T2D[length(brks_T2D)] <- 1.1 * brks_T2D[length(brks_T2D)]
dec_T2D <- cut(
  MY_T2D_WHITE$prob, 
  breaks = brks_T2D,
  include_lowest = TRUE
)

dec_T2D <- data.frame(y = MY_T2D_WHITE$M, pred = MY_T2D_WHITE$prob, dec = dec_T2D) %>%
  group_by(dec) %>%
  mutate(prob_obs = sum(y) / n(), 
         obs = sum(y),
         n_group = n(),
         mnpred = mean(pred),
         lower = lapply(sum(y), prop.test, n = n()), 
         upper = sapply(lower, function(x) x$conf.int[2]), 
         lower = sapply(lower, function(x) x$conf.int[1]))

## plot 
ggplot(dec_T2D, aes(x = mnpred, y = prob_obs)) +
  geom_point() +
  xlab("Mean predicted probability in each decile") +
  ylab("Proportion of MODY in each decile") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  #ylim(c(0, 1)) + xlim(c(0, 1)) +
  geom_errorbar(aes(ymin = lower, ymax = upper))

##Thresholds -----------------------------------------------------------------------------
### 5%
MY_T2D_WHITE_5PERC <- MY_T2D %>%
  summarise(totalover = sum(prob >= 0.05),
            ncasespickedup = sum(prob >= 0.05 & M ==1),
            PPV = (sum(prob >= 0.05 & M ==1)/sum(prob >= 0.05))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.05 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.05 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.05 & M == 0)/sum(prob < 0.05))*100,
            Sensitivity = (sum(prob >= 0.05 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.05 & M == 0)/sum(M==0))*100)
### 10%
MY_T2D_WHITE_10PERC <- MY_T2D %>%
  summarise(totalover = sum(prob >= 0.1),
            ncasespickedup = sum(prob >= 0.1 & M ==1),
            PPV = (sum(prob >= 0.1 & M ==1)/sum(prob >= 0.1))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.1 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.1 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.1 & M == 0)/sum(prob < 0.1))*100,
            Sensitivity = (sum(prob >= 0.1 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.1 & M == 0)/sum(M==0))*100)

### 20%
MY_T2D_WHITE_20PERC <- MY_T2D %>%
  summarise(totalover = sum(prob >= 0.2),
            ncasespickedup = sum(prob >= 0.2 & M ==1),
            PPV = (sum(prob >= 0.2 & M ==1)/sum(prob >= 0.2))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.2 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.2 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.2 & M == 0)/sum(prob < 0.2))*100,
            Sensitivity = (sum(prob >= 0.2 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.2 & M == 0)/sum(M==0))*100)

MY_T2D_WHITE_THRESHOLDS <- rbind(MY_T2D_WHITE_5PERC, MY_T2D_WHITE_10PERC, MY_T2D_WHITE_20PERC)
write_xlsx(MY_T2D_WHITE_THRESHOLDS,"MY_T2D_WHITE_THRESHOLDS_table.xlsx")