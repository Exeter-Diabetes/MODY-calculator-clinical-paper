########################################################################################################
#var_characteristics (variable characteristics)

#Function to produce characteristics tables
#On either a whole dataset
#Or by a grouping categorical variable
#With choice of how to summarise numeric variables 
#And summarise variable missingness

##########################################################################################
#Arguments:
# varlist (variable list) - string list of numeric variable names in dataset 
#                           (must be in interger or numeric class)
# varlist_cat (variable list categorical) - string list of categorical variable 
#                         names in dataset (must be character or factor class)
# dataset - name of dataset where in variables in varlist and varlist_cat are found
# numeric_option - choice of three ways of summerising numeric variables: 
#                     mean with 95% CI (meanCI), 
#                     mean with standard deviation (meanSD), 
#                     median with interquartile range (medianIQR)
# missingness - Whether want summary information of variable missingness 
#               (yes/TRUE is default)
# group - Grouping variable want to split data by and get summary data for 
#         accordingly
# table_name - name to save table to envrionment and excell as
#             if no user defined table name use default:
#             (name of dataset)_(name of grouping variable if provided)_summary_table
#             dataset_name_grouping_var_name_summary_table
# p_value_testing - whether want to add p-value/statistical testing to table
# stat_test_cont - statistical testing option for numeric variables
# stat_test_cat - statistical testing option for categorical variables

var_characteristics <- function(
    varlist = varlist, 
    varlist_cat = varlist_cat, 
    dataset = dataset, 
    numeric_option = c("meanCI", "meanSD", "medianIQR"), 
    missingness = TRUE, 
    group = NULL, 
    table_name = NULL,
    p_value_testing = TRUE,
    stat_test_cont = c("Wilcox"),
    stat_test_cat = c("Chi_square")
    ){
  
  #browser()
  
  #Argument checks 
  ##Checks for "dataset"
  if(missing(dataset) | is.null(dataset)) {stop("'dataset' needs to be provided.")}
  ##Checks for "varlist"
  if(missing(varlist) | is.null(varlist)) {stop("'varlist' (list of numeric variable names in dataset) needs to be provided.")}
  ##Checks for "varlist_cat"
  if(missing(varlist_cat) | is.null(varlist_cat)) {stop("'varlist_cat' (list of categorical variable names in dataset) needs to be provided.")}
  ##Checks for "numeric_option"
  if(missing(numeric_option) | is.null(numeric_option)) {stop("'numeric_option' needs to be provided.")}
  if(!(numeric_option %in% c("meanCI", "meanSD", "medianIQR"))) {stop("'numeric_option' needs to be: 'meanCI', 'meanSD', or 'medianIQR'")}
  ##Checks for "missingness"
  if(!(missingness %in% c(TRUE, FALSE))) {stop("'missingness' needs to be: TRUE / FALSE")}
  ##Checks for p_value_testing
  if(!(p_value_testing %in% c(TRUE, FALSE))) {stop("'p_value_testing' needs to be: TRUE / FALSE")}
  ##Checks for Statistical tests
  if(p_value_testing==TRUE & (missing(stat_test_cont) | is.null(stat_test_cont))) {stop("'stat_test_cont' needs to be provided")}
  if(p_value_testing==TRUE & (missing(stat_test_cat) | is.null(stat_test_cat))) {stop("'stat_test_cat' needs to be provided")}
  if(!missing(stat_test_cont) & !(stat_test_cont %in% c("Wilcox"))) {stop("'stat_test_cont' needs to be: 'Wilcox'")}
  if(!missing(stat_test_cat) & !(stat_test_cat %in% c("Chi_square"))) {stop("'stat_test_cat' needs to be: 'Chi_square'")}
  
  
  #Load libraries
  library(tidyverse)
  library(writexl)
  library(gdata)
  
  #save names as vectors
  dataset_name <- deparse(substitute(dataset))
  grouping_var_name <- deparse(substitute(group))
  name <- table_name
  num_option_name <- ifelse(numeric_option == "meanCI", "Mean (95% CI)",
                            ifelse(numeric_option == "meanSD", "Mean (SD)",
                                   "Median [IQR]"))
  #load external functions
  if(p_value_testing==TRUE){
    source("Functions/wilcox_table.R")
    source("Functions/chi_table.R")
  }
  
  source("Functions/rearrange_columns.R")
  #load internal functions
  

  #If no grouping variable provided
  if(is.null(group)){ 
    ##For categorical variables
    summaryTable_ncount <- data.frame(variable = NA, N = NA)
    summaryTable_ncount$variable <- "N"
    summaryTable_ncount$N <- nrow(dataset)
    for(var_cat in varlist_cat){
      sumtbl <- data.frame(variable = NA, category = NA, n_count = NA, perc = NA)  
      sumtbl$variable <- var_cat
      categories <- unique(dataset[[var_cat]])
      for (cat in categories){
        sumtbl$category <- as.character(cat)
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
    #colnames(summaryTable_ncount)[colnames(summaryTable_ncount) == "n_count"] <- "n"
    summaryTable_ncount <- as.data.frame(summaryTable_ncount)
    
    ##For numerical variables
    ###If numeric option is meanCI (mean & 95% Confidence Interval)  
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
      colnames(summaryTable_mCI)[colnames(summaryTable_mCI) == "mean_ci"] <- num_option_name
      
      summaryTable_count_num <- full_join(summaryTable_ncount, summaryTable_mCI)
    }
    
    ###If numeric option is meanSD (mean & standard deviation of the mean)  
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
      colnames(summaryTable_mSD)[colnames(summaryTable_mSD) == "mean_sd"] <- num_option_name
      summaryTable_mSD <- as.data.frame(summaryTable_mSD)
      summaryTable_count_num <- full_join(summaryTable_ncount, summaryTable_mSD)
    }
    
    ###If numeric option is medianIQR (median & interquartile range) 
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
      colnames(summaryTable_mIQR)[colnames(summaryTable_mIQR) == "median_IQR"] <- num_option_name
      summaryTable_mIQR <- as.data.frame(summaryTable_mIQR)
      summaryTable_count_num <- full_join(summaryTable_ncount, summaryTable_mIQR)
    }
    
    ##If adding missingness summary information
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
      colnames(summaryTable_nmiss)[colnames(summaryTable_nmiss) == "n_miss"] <- "n missing(%)"
      summaryTable_nmiss <- as.data.frame(summaryTable_nmiss)
      
      
      summaryTable_count_num_miss <- full_join(summaryTable_count_num, summaryTable_nmiss)
      #Save summary table to home environment
      #If user-defined summary table name not provided
      #use default (dataset_summary_table)
      if(is.null(table_name)) {
        #replace empty missing category rows with "Missing" string
        summaryTable_count_num_miss <-  summaryTable_count_num_miss %>%
          mutate(category = ifelse(summaryTable_count_num[[4]] == "NA (NA%)",
                                   NA,
                                   ifelse(!is.na(summaryTable_count_num[[4]]) & 
                                            is.na(category), 
                                          "Missing", 
                                          category)
          ),
          n_count = ifelse(n_count == "NA (NA%)",
                           N, n_count)) %>%
          select(-2)
        #for each category in group
        #make a new column that combines the category n column
        #and the numeric variable summary column
        #into one
        summaryTable_count_num_miss <-  summaryTable_count_num_miss %>%
            mutate(
              !!paste("n/", num_option_name, sep = " ") := 
                coalesce(!!sym("n_count"),
                         !!sym(num_option_name))
            ) %>%
            # Remove original 'n' columns for the group
            select(-3, -4)
        #reorder dataframe columns to be by group
        
        # Apply the function to rearrange the columns
        summaryTable_count_num_miss <- rearrange_columns(summaryTable_count_num_miss)
        
        #Save to global environment
        #and as Excell table
        #write_xlsx(summaryTable_count_num_miss, paste0(name,".xlsx"))
        summaryTable_count_num_miss <<- as.data.frame(summaryTable_count_num_miss)
        mv(from= "summaryTable_count_num_miss", 
           to = paste0(dataset_name, "_summary_table"), 
           envir = globalenv())
        write_xlsx(summaryTable_count_num_miss, 
                   paste0(dataset_name, "_summary_table.xlsx"))
        
      } else {
        #If user-defined table name provided
        #replace empty missing category rows with "Missing" string
        summaryTable_count_num_miss <-  summaryTable_count_num_miss %>%
          mutate(category = ifelse(summaryTable_count_num[[4]] == "NA (NA%)",
                                   NA,
                                   ifelse(!is.na(summaryTable_count_num[[4]]) & 
                                            is.na(category), 
                                          "Missing", 
                                          category)
          ),
          n_count = ifelse(n_count == "NA (NA%)",
                           N, n_count)) %>%
          select(-2)
        summaryTable_count_num_miss <- summaryTable_count_num_miss %>%
            mutate(
              !!paste("n/", num_option_name, sep = " ") := coalesce(!!sym("n_count"), 
                                                                       !!sym(num_option_name))
            ) %>%
            select(-3, -4)
        #reorder dataframe columns to be by group
        # Apply the function to rearrange the columns
        summaryTable_count_num_miss <- rearrange_columns(summaryTable_count_num_miss)
        
        #Save to global environment and Excell table 
        write_xlsx( summaryTable_count_num_miss, paste0(name,".xlsx"))
        summaryTable_count_num_miss <<- as.data.frame(summaryTable_count_num_miss)
        mv(from= "summaryTable_count_num_miss", to = paste0(name), envir = globalenv())
        #write_xlsx(summaryTable_GROUP_missing, paste0(name,".xlsx"))
      }
      } else {
      #If missingness variable false
      #Save summary table to home environment
      if(is.null(table_name)) {
        #replace empty missing category rows with "Missing" string
        summaryTable_count_num <-  summaryTable_count_num %>%
          mutate(category = ifelse(summaryTable_count_num[[4]] == "NA (NA%)",
                                   NA,
                                   ifelse(!is.na(summaryTable_count_num[[4]]) & 
                                            is.na(category), 
                                          "Missing", 
                                          category)
          ),
          n_count = ifelse(n_count == "NA (NA%)",
                           N, n_count)) %>%
          select(-2)
        #for each category in group
        #make a new column that combines the category n column
        #and the numeric variable summary column
        #into one
        summaryTable_count_num <- summaryTable_count_num %>%
            mutate(
              !!paste("n/", num_option_name, sep = " ") := 
                coalesce(!!sym(paste("n_count", sep = " ")),
                         !!sym(paste(num_option_name, sep = " ")))
            ) %>%
            # Remove original 'n' columns for the group
            select(-3, -4)
        #
        
        #Save to global environment
        #and as Excell table
        summaryTable_count_num <<- as.data.frame(summaryTable_count_num)
        mv(from= "summaryTable_count_num", to = paste0(dataset_name, "_summary_table"), envir = globalenv())
        write_xlsx(summaryTable_count_num, paste0(dataset_name, "_summary_table.xlsx"))
      } else {
        #If user-defined table name provided
        #summaryTable_GROUP <- full_join(summaryTable_GROUP, summaryTable_nmiss_group)
        #replace empty missing category rows with "Missing" string
        summaryTable_count_num <-  summaryTable_count_num %>%
          mutate(category = ifelse(summaryTable_count_num[[4]] == "NA (NA%)",
                                   NA,
                                   ifelse(!is.na(summaryTable_count_num[[4]]) & 
                                            is.na(category), 
                                          "Missing", 
                                          category)
                                   ),
                 n_count = ifelse(n_count == "NA (NA%)",
                                  N, n_count)) %>%
          select(-2)
          summaryTable_count_num <- summaryTable_count_num %>%
            mutate(
              !!paste("n/", num_option_name, sep = " ") := coalesce(!!sym(paste("n_count", sep = " ")), 
                                                                       !!sym(paste(num_option_name, sep = " ")))
            ) %>%
            select(-3, -4)
          #
          
          #
        summaryTable_count_num <<- as.data.frame(summaryTable_count_num)
        write_xlsx(summaryTable_count_num, paste0(name,".xlsx"))
        mv(from= "summaryTable_count_num", to = paste0(name), envir = globalenv())
      }
    }
  }
  
  #If grouping variable provided
  else{
    ##For categorical variables
    #browser()
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
          sumtbl$category <- as.character(cat)
          sumtbl$n <- sum(dataset1[[var]] == cat, na.rm = TRUE)
          sumtbl$perc <- (sum(dataset1[[var]] == cat, na.rm = TRUE)/length(dataset1[[var]]))*100
          #sumtbl$total <- nrow(dataset1[dataset1[[var]],])
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
      #summaryTable_ncount1$total <- NULL
      summaryTable_ncount_group <- left_join(summaryTable_ncount_group, summaryTable_ncount1)
    }
    #Save summary table to home environment
    summaryTable_ncount_group <- as.data.frame(summaryTable_ncount_group)
    
    if(p_value_testing == TRUE){
      if(stat_test_cat == "Chi_square"){
        chi_table(varlist_cat, dataset, outcome = group)
        
        chiTable <- chiTable %>%
          mutate(investigate = ifelse(p_value <= 0.05, 
                                      "Investigate", 
                                      "No difference between groups"))
      }
      if(stat_test_cont == "Wilcox"){
        wilcox_table(varlist, dataset, outcome = group)
        
        WilTable <- WilTable %>%
          mutate(investigate = ifelse(p_value <= 0.05, 
                                      "Investigate", 
                                      "No difference between groups"))
      }
    }
    
    ##For numeric variables
    ###If numeric option is meanCI (mean & 95% Confidence Interval) 
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
        colnames(summaryTable_mCI)[colnames(summaryTable_mCI) == "mean_ci"] <- paste(g, num_option_name)
        summaryTable_mCI_group <- left_join(summaryTable_mCI_group, summaryTable_mCI)  
      }
      #Save summary table to home environment
      summaryTable_GROUP <- full_join(summaryTable_ncount_group, summaryTable_mCI_group)
    }
    ###If numeric option is meanSD (mean & standard deviation)
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
        colnames(summaryTable_mSD)[colnames(summaryTable_mSD) == "mean_sd"] <- paste(g, num_option_name)
        summaryTable_mSD_group <- left_join(summaryTable_mSD_group, summaryTable_mSD)
      }
      #Save summary table to home environment
      summaryTable_GROUP <- full_join(summaryTable_ncount_group, summaryTable_mSD_group)
    }
    
    ###If numeric option is medianIQR (median & interquartile range)
    if(numeric_option == "medianIQR"){
      #browser()
      summaryTable_mIQR_group <- data.frame(variable = varlist)
      for(g in groups){
        summaryTable_mIQR <- data.frame()
        dataset1 <- dataset[dataset[[group]]== g & !is.na(dataset[[group]]),]
        for(var in varlist){
          sumtbl <- data.frame(variable = NA, median = NA, Q1 = NA, Q3 = NA)  
          sumtbl$variable <- var
          print(var)
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
        colnames(summaryTable_mIQR)[colnames(summaryTable_mIQR) == "median_IQR"] <- paste(g, num_option_name)
        summaryTable_mIQR_group <- left_join(summaryTable_mIQR_group, summaryTable_mIQR)
      }
      #Save summary table to home environment
      summaryTable_GROUP <- full_join(summaryTable_ncount_group, summaryTable_mIQR_group)
    }
    
    ##If adding missingness summary information
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
        summaryTable_nmiss_group <- left_join(summaryTable_nmiss_group, 
                                              summaryTable_nmiss)
      }
      #Save summary table to home environment
      #if no user defined table name use default 
      #dataset_name_grouping_var_name_summary_table
      if(is.null(table_name)) {
        summaryTable_GROUP_missing <- full_join(summaryTable_GROUP, 
                                                summaryTable_nmiss_group)
        #replace empty missing category rows with "Missing" string
        summaryTable_GROUP_missing <- summaryTable_GROUP_missing %>%
          mutate(category = ifelse(!is.na(summaryTable_GROUP_missing[[3]]) & 
                                     is.na(category), "Missing", category))
        #for each category in group
        #make a new column that combines the category n column
        #and the numeric variable summary column
        #into one
        groups <- unique(dataset[[group]])
        groups <- groups[ !groups == 'NA']
        for(g in groups) {
          size_pattern <- paste(c(g, num_option_name), collapse = "|")
          summaryTable_GROUP_missing <- summaryTable_GROUP_missing %>%
            mutate(
              !!paste(g, "n/", num_option_name, sep = " ") := 
                coalesce(!!sym(paste(g, "n", sep = " ")),
                         !!sym(paste(g, num_option_name, sep = " ")))
            ) %>%
            # Remove original 'n' columns for the group
            select(-matches(paste0("^", g, "( n| (", num_option_name, "))$", 
                                   sep = ""))) 
          
        }
        #Remove the numeric variable summary columns
        summaryTable_GROUP_missing <<- summaryTable_GROUP_missing %>%
          select(-3, -4)
        print(colnames(summaryTable_GROUP_missing))
        #reorder dataframe columns to be by group
        # Apply the function to rearrange the columns
        summaryTable_GROUP_missing <- rearrange_columns(summaryTable_GROUP_missing, groups, num_option_name)
        #
        if(p_value_testing == TRUE){
        table_pvalues <- rbind(chiTable, WilTable)
        summaryTable_GROUP_missing <- left_join(summaryTable_GROUP_missing, table_pvalues)
      } 
        #Save to global environment
        #and as Excell table
        write_xlsx(summaryTable_GROUP_missing, paste0(name,".xlsx"))
        summaryTable_GROUP_missing <<- as.data.frame(summaryTable_GROUP_missing)
        mv(from= "summaryTable_GROUP_missing", 
           to = paste0(dataset_name, "_", 
                       grouping_var_name, "_summary_table"), 
           envir = globalenv())
        write_xlsx(summaryTable_GROUP_missing, 
                   paste0(dataset_name,"_",
                          grouping_var_name,"_summary_table.xlsx"))
        
      } else {
        #If user-defined table name provided
        summaryTable_GROUP_missing <- full_join(summaryTable_GROUP, summaryTable_nmiss_group)
        #replace empty missing category rows with "Missing" string
        summaryTable_GROUP_missing <- summaryTable_GROUP_missing %>%
          mutate(category = ifelse(!is.na(summaryTable_GROUP_missing[[3]]) & is.na(category), "Missing", category))
        groups <- unique(dataset[[group]])
        groups <- groups[ !groups == 'NA']
        for(g in groups) {
          size_pattern <- paste(c(g, num_option_name), collapse = "|")
          summaryTable_GROUP_missing <- summaryTable_GROUP_missing %>%
            mutate(
              !!paste(g, "n/", num_option_name, sep = " ") := coalesce(!!sym(paste(g, "n", sep = " ")), 
                                                                       !!sym(paste(g, num_option_name, sep = " ")))
            ) #%>%
            #select(-matches(paste0("^", g, "( n| (", num_option_name, "))$", sep = ""))) 
           # Remove original 'n' and 'size' columns for the group
        }
        #reorder dataframe columns to be by group
        #print(colnames(summaryTable_GROUP_missing))
        # Apply the function to rearrange the columns
        summaryTable_GROUP_missing <- rearrange_columns(summaryTable_GROUP_missing, groups, num_option_name)
        #
        if(p_value_testing == TRUE){
          table_pvalues <- rbind(chiTable, WilTable)
          summaryTable_GROUP_missing <- left_join(summaryTable_GROUP_missing, table_pvalues)
        } 
        #Save to global environment and Excell table 
        write_xlsx(summaryTable_GROUP_missing, paste0(name,".xlsx"))
        summaryTable_GROUP_missing <<- as.data.frame(summaryTable_GROUP_missing)
        mv(from= "summaryTable_GROUP_missing", to = paste0(name), envir = globalenv())
        #write_xlsx(summaryTable_GROUP_missing, paste0(name,".xlsx"))
      }
    } else {
      #Save summary table to home environment
      if(is.null(table_name)) {
        #replace empty missing category rows with "Missing" string
        summaryTable_GROUP <- summaryTable_GROUP %>%
          mutate(category = ifelse(!is.na(summaryTable_GROUP[[3]]) & 
                                     is.na(category), "Missing", category))
        #for each category in group
        #make a new column that combines the category n column
        #and the numeric variable summary column
        #into one
        groups <- unique(dataset[[group]])
        groups <- groups[ !groups == 'NA']
        for(g in groups) {
          size_pattern <- paste(c(g, num_option_name), collapse = "|")
          summaryTable_GROUP <- summaryTable_GROUP %>%
            mutate(
              !!paste(g, "n/", num_option_name, sep = " ") := 
                coalesce(!!sym(paste(g, "n", sep = " ")),
                         !!sym(paste(g, num_option_name, sep = " ")))
            ) %>%
            # Remove original 'n' columns for the group
            select(-matches(paste0("^", g, "( n| (", num_option_name, "))$", 
                                   sep = ""))) 
        }
        #change this to subsetting in rearrange_columns
        summaryTable_GROUP <- summaryTable_GROUP %>%
          select(-3, -4)
        #
        if(p_value_testing == TRUE){
          table_pvalues <- rbind(chiTable, WilTable_SR)
          summaryTable_GROUP <- left_join(summaryTable_GROUP, table_pvalues)
        } 
        #Save to global environment
        #and as Excell table
        summaryTable_GROUP <<- as.data.frame(summaryTable_GROUP)
        mv(from= "summaryTable_GROUP", to = paste0(dataset_name, "_", grouping_var_name, "_summary_table"), envir = globalenv())
        write_xlsx(summaryTable_GROUP, paste0(dataset_name,"_",grouping_var_name,"_summary_table.xlsx"))
      } else {
        #If user-defined table name provided
        #summaryTable_GROUP <- full_join(summaryTable_GROUP, summaryTable_nmiss_group)
        #replace empty missing category rows with "Missing" string
        summaryTable_GROUP <- summaryTable_GROUP %>%
          mutate(category = ifelse(!is.na(summaryTable_GROUP[[3]]) & is.na(category), "Missing", category))
        groups <- unique(dataset[[group]])
        groups <- groups[ !groups == 'NA']
        for(g in groups) {
          size_pattern <- paste(c(g, num_option_name), collapse = "|")
          summaryTable_GROUP <- summaryTable_GROUP %>%
            mutate(
              !!paste(g, "n/", num_option_name, sep = " ") := coalesce(!!sym(paste(g, "n", sep = " ")), 
                                                                       !!sym(paste(g, num_option_name, sep = " ")))
            ) %>%
            select(-matches(paste0("^", g, "( n| (", num_option_name, "))$", sep = ""))) 
          # Remove original 'n' and 'size' columns for the group
        }
        summaryTable_GROUP <- summaryTable_GROUP %>%
          select(-3, -4)
        #
        if(p_value_testing == TRUE){
          table_pvalues <- rbind(chiTable, WilTable)
          summaryTable_GROUP <- left_join(summaryTable_GROUP, table_pvalues)
        } 
        #
        summaryTable_GROUP <<- as.data.frame(summaryTable_GROUP)
        write_xlsx(summaryTable_GROUP, paste0(name,".xlsx"))
        mv(from= "summaryTable_GROUP", to = paste0(name), envir = globalenv())
      }
    }
  }
}





