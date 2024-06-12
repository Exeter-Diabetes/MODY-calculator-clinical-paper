## FUNCTION WITHOUT GROUPBY OPTION
var_characteristics <- function(varlist = varlist, varlist_cat = varlist_cat, dataset = dataset, numeric_option = c("meanCI", "meanSD", "medianIQR"), missingness = TRUE, group = NULL){
  #browser()
  
  if(is.null(group)){ 
    #FOR CATEGORICAL VARIABLES
    summaryTable_ncount <- data.frame(variable = NA, N = NA)
    summaryTable_ncount$variable <- "N"
    summaryTable_ncount$N <- nrow(dataset)
    for(var_cat in varlist_cat){
      sumtbl <- data.frame(variable = NA, category = NA, n_count = NA, perc = NA)  
      sumtbl$variable <- var_cat
      categories <- unique(dataset[[var_cat]])
      for (cat in categories){
        
        sumtbl$category <- cat
        sumtbl$n <- sum(dataset[[var_cat]] == cat, na.rm = TRUE)
        sumtbl$perc <- (sum(dataset[[var_cat]] == cat, na.rm = TRUE)/length(dataset[[var_cat]]))*100
        sumtbl$total <- nrow(dataset[dataset[[var_cat]],])
        summaryTable_ncount <- bind_rows(summaryTable_ncount, sumtbl)
      }
      
    }
    #combine mean, lci and uci columns into 1 (mean_ci)
    summaryTable_ncount$perc_1 <- paste0("(", round(summaryTable_ncount$perc, 2),"%",")")  
    summaryTable_ncount$n_count <- paste0(summaryTable_ncount$n, " ", summaryTable_ncount$perc_1)
    summaryTable_ncount$perc <- NULL 
    summaryTable_ncount$n <- NULL
    summaryTable_ncount$perc_1 <- NULL
    summaryTable_ncount$total <- NULL
    
    summaryTable_ncount <- as.data.frame(summaryTable_ncount)
    
    #IF WANT MEAN & 95% CI FOR NUMERIC VARIABLES IN VARLIST  
    if(numeric_option == "meanCI"){
      summaryTable_mCI <- data.frame()
      for(var in varlist){
        
        sumtbl <- data.frame(variable = NA, mean_ci = NA)  
        sumtbl$variable <- var
        sumtbl$mean_ci <- list(mean_cl_normal(dataset[[var]]) %>%
                                 rename(mean = y, lci = ymin, uci = ymax)) 
        summaryTable_mCI <- rbind(summaryTable_mCI, sumtbl)
      }
      summaryTable_mCI <- as.data.frame(summaryTable_mCI) %>%
        unnest(cols = mean_ci) 
      #convert to two decimals
      summaryTable_mCI <- summaryTable_mCI %>%
        mutate_if(is.numeric, round, digits = 2)
      #combine mean, lci and uci columns into 1 (mean_ci)
      summaryTable_mCI$uci <- paste(summaryTable_mCI$uci, ")")  
      summaryTable_mCI$lci <- paste("(", summaryTable_mCI$lci) 
      summaryTable_mCI <- summaryTable_mCI %>%
        unite(lci_uci, c(lci, uci), sep = " ,", remove = TRUE) %>%
        unite(mean_ci, c(mean,lci_uci), sep = " ", remove = TRUE)
      
      summaryTable_count_num <- full_join(summaryTable_ncount, summaryTable_mCI)
      
      
    }
    #IF WANT MEAN & SD FOR NUMERIC VARIABLES IN VARLIST  
    if(numeric_option == "meanSD"){
      summaryTable_mSD <- data.frame()
      for(var in varlist){
        
        sumtbl <- data.frame(variable = NA, mean = NA, sd = NA)  
        sumtbl$variable <- var
        sumtbl$mean <-  mean(dataset[[var]], na.rm = TRUE)
        sumtbl$sd <- sd(dataset[[var]], na.rm = TRUE)
        summaryTable_mSD <- rbind(summaryTable_mSD, sumtbl)
      }
      
      #convert to two decimals
      summaryTable_mSD <- summaryTable_mSD %>%
        mutate_if(is.numeric, round, digits = 2)
      #combine mean, lci and uci columns into 1 (mean_ci)
      summaryTable_mSD$sd <- paste0("(", round(summaryTable_mSD$sd, digits = 2), ")")  
      summaryTable_mSD$mean_sd <- paste0(round(summaryTable_mSD$mean, digits = 2), " ", summaryTable_mSD$sd)
      summaryTable_mSD$mean <- NULL 
      summaryTable_mSD$sd <- NULL
      summaryTable_mSD <- as.data.frame(summaryTable_mSD)
      summaryTable_count_num <- full_join(summaryTable_ncount, summaryTable_mSD)
    }
    #IF WANT MEDIAN AND IQR FOR NUMERIC VARIABLES IN VARLIST  
    if(numeric_option == "medianIQR"){
      summaryTable_mIQR <- data.frame()
      for(var in varlist){
        sumtbl <- data.frame(variable = NA, median = NA, Q1 = NA, Q3 = NA)  
        sumtbl$variable <- var
        sumtbl$median <-  median(dataset[[var]], na.rm = TRUE)
        sumtbl$Q1 <- quantile(dataset[[var]], na.rm = TRUE, 0.25)
        sumtbl$Q3 <- quantile(dataset[[var]], na.rm = TRUE, 0.75)
        summaryTable_mIQR <- rbind(summaryTable_mIQR, sumtbl)
      }
      
      #convert to two decimals
      summaryTable_mIQR <- summaryTable_mIQR %>%
        mutate_if(is.numeric, round, digits = 2)
      #combine mean, lci and uci columns into 1 (mean_ci)
      summaryTable_mIQR$IQR <- paste0("[", round(summaryTable_mIQR$Q1, digits = 2),",", round(summaryTable_mIQR$Q3, 2), "]")  
      summaryTable_mIQR$median_IQR <- paste0(round(summaryTable_mIQR$median, digits = 2), " ", summaryTable_mIQR$IQR)
      summaryTable_mIQR$median <- NULL 
      summaryTable_mIQR$Q1 <- NULL
      summaryTable_mIQR$Q3 <- NULL
      summaryTable_mIQR$IQR <- NULL
      summaryTable_mIQR <- as.data.frame(summaryTable_mIQR)
      summaryTable_count_num <- full_join(summaryTable_ncount, summaryTable_mIQR)
    }
    #IF WANT TO INCLUDE MISSINGNESS INFORMATION
    if(missingness == TRUE){
      summaryTable_nmiss <- data.frame()
      varlist_all <- c(varlist, varlist_cat)
      for(var_all in varlist_all){
        
        sumtbl <- data.frame(variable = NA, n_miss = NA, perc = NA)  
        sumtbl$variable <- var_all
        
        sumtbl$n <- sum(is.na(dataset[[var_all]]) == TRUE, na.rm = TRUE)
        #sumtbl$perc <- (sum(is.na(dataset[[var_all]]) == TRUE, na.rm = TRUE)/nrow(dataset[[var_all]]))*100
        sumtbl$perc <- (sum(is.na(dataset[[var_all]]) == TRUE, na.rm = TRUE)/length(dataset[[var_all]]))*100
        
        summaryTable_nmiss <- rbind(summaryTable_nmiss, sumtbl)
      }
      
      #combine mean, lci and uci columns into 1 (mean_ci)
      summaryTable_nmiss$perc_1 <- paste0("(", round(summaryTable_nmiss$perc, 2),"%",")")  
      summaryTable_nmiss$n_miss <- paste0(summaryTable_nmiss$n, " ", summaryTable_nmiss$perc_1)
      summaryTable_nmiss$perc <- NULL 
      summaryTable_nmiss$n <- NULL
      summaryTable_nmiss$perc_1 <- NULL
      summaryTable_nmiss <- as.data.frame(summaryTable_nmiss)
      
      summaryTable_count_num_miss <<- full_join(summaryTable_count_num, summaryTable_nmiss)
    }
  }
  else{
    #FOR CATEGORICAL VARIABLES
    summaryTable_ncount_group <- data.frame(variable = varlist_cat)
    groups <- unique(dataset[[group]])
    groups <- groups[ !groups == 'NA']
    for (g in groups) {
      summaryTable_ncount <- data.frame()
      dataset1 <- dataset[dataset[[group]]== g & !is.na(dataset[[group]]),]
      for(var in varlist_cat){
        sumtbl <- data.frame(variable = NA, category = NA, n = NA, perc = NA)  
        sumtbl$variable <- var
        categories <- unique(dataset[[var]])
        
        for (cat in categories){
          sumtbl$category <- cat
          sumtbl$n <- sum(dataset1[[var]] == cat, na.rm = TRUE)
          sumtbl$perc <- (sum(dataset1[[var]] == cat, na.rm = TRUE)/length(dataset1[[var]]))*100
          sumtbl$total <- nrow(dataset1[dataset1[[var]],])
          summaryTable_ncount <- rbind(summaryTable_ncount, sumtbl)
        }
        summaryTable_ncount1 <- left_join(summaryTable_ncount_group, summaryTable_ncount)
      }
      
      summaryTable_ncount1$perc_1 <- paste0("(", round(summaryTable_ncount1$perc, 2),"%",")")  
      summaryTable_ncount1$n_perc <- paste0(summaryTable_ncount1$n, " ", summaryTable_ncount1$perc_1)
      colnames(summaryTable_ncount1)[colnames(summaryTable_ncount1) == "n_perc"] <- paste(g, "n")
      summaryTable_ncount1$perc <- NULL 
      summaryTable_ncount1$n <- NULL
      summaryTable_ncount1$perc_1 <- NULL
      summaryTable_ncount1$total <- NULL
      summaryTable_ncount_group <- left_join(summaryTable_ncount_group, summaryTable_ncount1)
    }
    summaryTable_ncount_group <<- as.data.frame(summaryTable_ncount_group)
    
    ## FOR NUMERIC VARIABLES
    #IF WANT MEAN & 95% CI 
    if(numeric_option == "meanCI"){
      summaryTable_mCI_group <- data.frame(variable = varlist)
      for(g in groups){
        summaryTable_mCI <- data.frame()
        dataset1 <- dataset[dataset[[group]]== g & !is.na(dataset[[group]]),]
        for(var in varlist){
          
          sumtbl <- data.frame(variable = NA, mean_ci = NA)  
          sumtbl$variable <- var
          sumtbl$mean_ci <- list(mean_cl_normal(dataset1[[var]]) %>%
                                   rename(mean = y, lci = ymin, uci = ymax)) 
          summaryTable_mCI <- rbind(summaryTable_mCI, sumtbl)
        }
        summaryTable_mCI <- as.data.frame(summaryTable_mCI) %>%
          unnest(cols = mean_ci) 
        #convert to two decimals
        summaryTable_mCI <- summaryTable_mCI %>%
          mutate_if(is.numeric, round, digits = 2)
        #combine mean, lci and uci columns into 1 (mean_ci)
        summaryTable_mCI$uci <- paste(summaryTable_mCI$uci, ")")  
        summaryTable_mCI$lci <- paste("(", summaryTable_mCI$lci) 
        summaryTable_mCI <- summaryTable_mCI %>%
          unite(lci_uci, c(lci, uci), sep = " ,", remove = TRUE) %>%
          unite(mean_ci, c(mean,lci_uci), sep = " ", remove = TRUE)
        colnames(summaryTable_mCI)[colnames(summaryTable_mCI) == "mean_ci"] <- paste(g, "mean_ci")
        summaryTable_mCI_group <- left_join(summaryTable_mCI_group, summaryTable_mCI)  
      }
      summaryTable_GROUP <<- full_join(summaryTable_ncount_group, summaryTable_mCI_group)
      
    }
    #IF WANT MEAN & SD
    if(numeric_option == "meanSD"){
      summaryTable_mSD_group <- data.frame(variable = varlist)
      for(g in groups){
        summaryTable_mSD <- data.frame()
        dataset1 <- dataset[dataset[[group]]== g & !is.na(dataset[[group]]),]
        for(var in varlist){
          
          sumtbl <- data.frame(variable = NA, mean = NA, sd = NA)  
          sumtbl$variable <- var
          sumtbl$mean <-  mean(dataset1[[var]], na.rm = TRUE)
          sumtbl$sd <- sd(dataset1[[var]], na.rm = TRUE)
          summaryTable_mSD <- rbind(summaryTable_mSD, sumtbl)
        }
        #convert to two decimals
        summaryTable_mSD <- summaryTable_mSD %>%
          mutate_if(is.numeric, round, digits = 2)
        #combine mean, lci and uci columns into 1 (mean_ci)
        summaryTable_mSD$sd <- paste0("(", round(summaryTable_mSD$sd, digits = 2), ")")  
        summaryTable_mSD$mean_sd <- paste0(round(summaryTable_mSD$mean, digits = 2), " ", summaryTable_mSD$sd)
        summaryTable_mSD$mean <- NULL 
        summaryTable_mSD$sd <- NULL
        colnames(summaryTable_mSD)[colnames(summaryTable_mSD) == "mean_sd"] <- paste(g, "mean_sd")
        summaryTable_mSD_group <- left_join(summaryTable_mSD_group, summaryTable_mSD)
      }
      
      summaryTable_GROUP <<- full_join(summaryTable_ncount_group, summaryTable_mSD_group)
    }
    
    #IF WANT MEDIAN & IQR
    browser()
    if(numeric_option == "medianIQR"){
      summaryTable_mIQR_group <- data.frame(variable = varlist)
      for(g in groups){
        summaryTable_mIQR <- data.frame()
        dataset1 <- dataset[dataset[[group]]== g & !is.na(dataset[[group]]),]
        for(var in varlist){
          sumtbl <- data.frame(variable = NA, median = NA, Q1 = NA, Q3 = NA)  
          sumtbl$variable <- var
          sumtbl$median <-  median(dataset1[[var]], na.rm = TRUE)
          sumtbl$Q1 <- quantile(dataset1[[var]], na.rm = TRUE, 0.25)
          sumtbl$Q3 <- quantile(dataset1[[var]], na.rm = TRUE, 0.75)
          summaryTable_mIQR <- rbind(summaryTable_mIQR, sumtbl)
        }
        
        #convert to two decimals
        summaryTable_mIQR <- summaryTable_mIQR %>%
          mutate_if(is.numeric, round, digits = 2)
        #combine mean, lci and uci columns into 1 (mean_ci)
        summaryTable_mIQR$IQR <- paste0("[", round(summaryTable_mIQR$Q1, digits = 2),",", round(summaryTable_mIQR$Q3, 2), "]")  
        summaryTable_mIQR$median_IQR <- paste0(round(summaryTable_mIQR$median, digits = 2), " ", summaryTable_mIQR$IQR)
        summaryTable_mIQR$median <- NULL 
        summaryTable_mIQR$Q1 <- NULL
        summaryTable_mIQR$Q3 <- NULL
        summaryTable_mIQR$IQR <- NULL
        colnames(summaryTable_mIQR)[colnames(summaryTable_mIQR) == "median_IQR"] <- paste(g, "median_IQR")
        summaryTable_mIQR_group <- left_join(summaryTable_mIQR_group, summaryTable_mIQR)
      }
      summaryTable_GROUP <<- full_join(summaryTable_ncount_group, summaryTable_mIQR_group)
    }
    
    ## IF WANT TO INCLUDE MISSINGNESS
    if(missingness == TRUE){
      varlist_all <- c(varlist, varlist_cat)
      summaryTable_nmiss_group <- data.frame(variable = varlist_all)
      for(g in groups){
        summaryTable_nmiss <- data.frame()
        dataset1 <- dataset[dataset[[group]]== g & !is.na(dataset[[group]]),]
        for(var_all in varlist_all){
          #if(var_all == "bmi") browser()
          sumtbl <- data.frame(variable = NA, n_miss = NA, perc = NA)  
          sumtbl$variable <- var_all
          
          sumtbl$n <- sum(is.na(dataset1[[var_all]]) == TRUE, na.rm = TRUE)
          #sumtbl$perc <- (sum(is.na(dataset[[var]]) == TRUE, na.rm = TRUE)/nrow(dataset[[var]]))*100
          sumtbl$perc <- (sum(is.na(dataset1[[var_all]]) == TRUE, na.rm = TRUE)/length(dataset1[[var_all]]))*100
          
          summaryTable_nmiss <- rbind(summaryTable_nmiss, sumtbl)
        }
        
        
        #combine mean, lci and uci columns into 1 (mean_ci)
        summaryTable_nmiss$perc_1 <- paste0("(", round(summaryTable_nmiss$perc, 2),"%",")")  
        summaryTable_nmiss$n_miss <- paste0(summaryTable_nmiss$n, " ", summaryTable_nmiss$perc_1)
        summaryTable_nmiss$perc <- NULL 
        summaryTable_nmiss$n <- NULL
        summaryTable_nmiss$perc_1 <- NULL
        colnames(summaryTable_nmiss)[colnames(summaryTable_nmiss) == "n_miss"] <- paste(g, "n missing(%)")
        summaryTable_nmiss_group <- left_join(summaryTable_nmiss_group, summaryTable_nmiss)
        
      }
      summaryTable_GROUP_missing <<- full_join(summaryTable_GROUP, summaryTable_nmiss_group)
    }
  }
}




#var_characteristics(varlist, varlist_cat, dataset, numeric_option = "medianIQR", group = "insulinRequire3")
#var_characteristics(varlist, varlist_cat, dataset, numeric_option = "medianIQR")
