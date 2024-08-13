#UNITED PAEDIATRIC VALIDATION ---------------------------------------------------------------
# load libraries ----------------------------------------------------------------------------
library(nimble)
library(rms)
library(tidyverse)
library(writexl)

## load functions ------------------------------------------------------------------------
# needed for generating data 
source("data/create_data.R")

#load in UNITED PAEDIATRIC datasets -----------------------------------------------------
#M = MODY status (NA = not tested for MODY (i.e. T=1), 1 = tested for MODY (T=0) and positive, 0 = tested for MODY (T=0) and negative)
#C = c-peptide status (1 = UCPCR >= 0.2; 0 = UCPCR < 0.2)
#A = antibody status (1 = 1+ positive antibody, 0 = all antibodies tested negative)
#pardm = parent with diabetes (1 = 1+ parent affected, 0 = no parents affected)
#insoroha = currently treated with insulin or oha (tablets) (1= currently treated with either, 0 = not treated with either)
#T = biomarker status (1 = cpeptide negative (UCPCR < 0.2) or antibody positive (A =1), 0 = cpeptide positive (UCPCR >=0.2) AND antibody negative (A=0))
dataset.UNITED_type1p <- create_data(dataset = "united t1d pediatrics", commonmody = FALSE)

#need to change M=NA to M=0
dataset.UNITED_type1p <- dataset.UNITED_type1p %>%
  mutate(M = ifelse(is.na(M) == TRUE, 0, M))
#checked if worked: should have M=1 (n=7) & M=0 (n=1164)
table(dataset.UNITED_type1p$M)


#load predictions --------------------------------------------------------------------------------
#load T1D
predictions_dataset.UNITED_type1_young_all_genes_with_T <- readRDS("model_predictions/predictions_dataset.UNITED_type1_young_all_genes_with_T.rds")
#JOIN PREDICTIONS AND DATA
UNITED_type1p <- cbind(dataset.UNITED_type1p, predictions_dataset.UNITED_type1_young_all_genes_with_T)


## Early-insulin-treated model analysis -----------------------------------------------------------------------
UNITED_type1p %>%
  group_by(M) %>%
  summarise(mean = mean(prob),
            meanLCI = mean(LCI),
            meanUCI = mean(UCI),
            n= n())


##BOXPLOTS ----------------------------------------------------------------------------------------
ggplot(UNITED_type1p, aes(x = factor(M, levels = c(0, 1), labels = c("Negative", "Positive")), y = prob)) +
  geom_boxplot() +
  labs(x = "MODY Status", y = "Probability") +
  theme_bw()

##ROC -----------------------------------------------------------------------------------------

roc_UNITED_type1p_model <- roc(UNITED_type1p$M, UNITED_type1p$prob, plot = TRUE, print.thres = "best", print.auc = TRUE)
model_UNITED_type1p_pr <- coords(roc_UNITED_type1p_model, x = "best", ret=c("threshold", "specificity", "sensitivity", "accuracy", "precision", "recall", "ppv", "npv"), transpose = FALSE)
#model_pr <- model_pr[c(-1,-3),]


sum_UNITED_type1ptable$threshold <- model_UNITED_type1p_pr$threshold
sum_UNITED_type1ptable$ROCAUC <- roc_UNITED_type1p_model$auc
sum_UNITED_type1ptable$accuracy <- model_UNITED_type1p_pr$accuracy
sum_UNITED_type1ptable$sensitivity <- model_UNITED_type1p_pr$sensitivity
sum_UNITED_type1ptable$specificity <- model_UNITED_type1p_pr$specificity
sum_UNITED_type1ptable$PPV <- model_UNITED_type1p_pr$ppv
sum_UNITED_type1ptable$NPV <- model_UNITED_type1p_pr$npv

print(roc_UNITED_type1p_model)

UNITED_type1p_ROC_table <- as.data.frame(sum_UNITED_type1ptable)
write_xlsx(UNITED_type1p_ROC_table,"UNITED_type1p_ROC_table.xlsx")

## Calibration ---------------------------------------------------------------------------------------
brks_t1d <- unique(quantile(UNITED_type1p$prob, prob = seq(0, 1, by = 0.2), na.rm = TRUE))
brks_t1d[1] <- 0.99 * brks_t1d[1]
brks_t1d[length(brks_t1d)] <- 1.1 * brks_t1d[length(brks_t1d)]
dec_t1d <- cut(
  UNITED_type1p$prob, 
  breaks = brks_t1d,
  include_lowest = TRUE
)

dec_t1d <- data.frame(y = UNITED_type1p$M, pred = UNITED_type1p$prob, dec = dec_t1d) %>%
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
UNITED_type1p_WHITE_5PERC <- UNITED_type1p %>%
  summarise(totalover = sum(prob >= 0.05),
            ncasespickedup = sum(prob >= 0.05 & M ==1),
            PPV = (sum(prob >= 0.05 & M ==1)/sum(prob >= 0.05))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.05 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.05 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.05 & M == 0)/sum(prob < 0.05))*100,
            Sensitivity = (sum(prob >= 0.05 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.05 & M == 0)/sum(M==0))*100)
### 10%
UNITED_type1p_WHITE_10PERC <- UNITED_type1p %>%
  summarise(totalover = sum(prob >= 0.1),
            ncasespickedup = sum(prob >= 0.1 & M ==1),
            PPV = (sum(prob >= 0.1 & M ==1)/sum(prob >= 0.1))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.1 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.1 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.1 & M == 0)/sum(prob < 0.1))*100,
            Sensitivity = (sum(prob >= 0.1 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.1 & M == 0)/sum(M==0))*100)

### 20%
UNITED_type1p_WHITE_20PERC <- UNITED_type1p %>%
  summarise(totalover = sum(prob >= 0.2),
            ncasespickedup = sum(prob >= 0.2 & M ==1),
            PPV = (sum(prob >= 0.2 & M ==1)/sum(prob >= 0.2))*100,
            nmissedcases = sum(M==1) - sum(prob >= 0.2 & M ==1),
            Missedcases = ((sum(M==1) - sum(prob >= 0.2 & M ==1))/sum(M==1))*100,
            NPV = (sum(prob < 0.2 & M == 0)/sum(prob < 0.2))*100,
            Sensitivity = (sum(prob >= 0.2 & M ==1)/sum(M==1))*100,
            Specificity = (sum(prob < 0.2 & M == 0)/sum(M==0))*100)

UNITED_type1p_WHITE_THRESHOLDS <- rbind(UNITED_type1p_WHITE_5PERC, UNITED_type1p_WHITE_10PERC, UNITED_type1p_WHITE_20PERC)
write_xlsx(UNITED_type1p_WHITE_THRESHOLDS,"UNITED_type1p_WHITE_THRESHOLDS_table.xlsx")
