###
# Funciton to make data spaces over range of values
###


make_prediction_data <- function(data, 
                                 col_names = NULL.
                                 len.out = 1000,
                                 constrain = T) {
  
  if(is.null(col_names)) {
    col_names = names(data)
  }
  
  data <- data |> 
    select(all_of(col_names))
  
  # check for dates
  if(any(sapply(data, is.Date))) {
    date_cols <- which(sapply(data, is.Date))
    data <- data[,-c(date_cols)]
    warning('Date columns are not compatible and were removed. Convert to numeric if interested.')
  }
  
  # check for posixct
  if(any(sapply(data, is.POSIXct))) {
    date_cols <- which(sapply(data, is.POSIXct))
    data <- data[,-c(date_cols)]
    warning('Date columns are not compatible and were removed. Convert to numeric if interested.')
  }
  
  ## Make grid of character & factor
  if(any(sapply(data, is.factor))) {
    fact_cols <- which(sapply(data, is.factor))
    data[fact_cols] <- lapply(data[fact_cols], as.character)
    warning('Factors are not currently supported, converted to char')
  }
  
  
  ##  For when character values are present
  char_cols <- which(sapply(data, is.character))
  if(length(char_cols) > 0) {
    hold_list <- list() # init list
    
    for(col in names(char_cols)) {
      hold_list[[col]] <- unique(data[[col]])
    }
    #sequence is important here
    udf <- expand.grid(hold_list) # make unique parings df
    unique_pairs <- nrow(outdf)
    # for len.out
    outdf <- udf[rep(1:nrow(udf), each = len.out),]
    row.names(outdf) <- NULL
    
    # simple approach 
    if(!constrain) {
      
    } else {
      
    }
    
  } else {
    # this is more simple, there are no character
  }

  
}