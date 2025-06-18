########################################################################################

######################################################################################

rearrange_columns <- function(df, groups = NULL, num_option_name = num_option_name) {
  #browser()
  if(is.null(groups)){
    new_order <- c("variable", "category")
    
    # Loop over each group to add columns in the "x" and "y" alternating pattern
    new_order <- c(new_order, paste("n/", num_option_name, sep = " "), paste("n missing(%)", sep = " "))
    
    # Reorder columns
    df <- df[, new_order]
    return(df)
  }else{
    # Start with the initial columns ("var" and "cat")
    new_order <- c("variable", "category")
    
    # Loop over each group to add columns in the "x" and "y" alternating pattern
    for(g in groups) {
      new_order <- c(new_order, paste(g, "n/", num_option_name, sep = " "), paste(g, "n missing(%)", sep = " "))
      #print(g)
      #print(new_order)
    }
    #print(new_order)
    df <- df %>%
      dplyr::select(all_of(new_order))
    #print(colnames(df))
    # Reorder columns
    df <- df[, new_order]
    return(df)
  }
}