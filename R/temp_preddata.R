###
# Funciton to make data spaces over range of values
###

#' Note that this is useful for generating all data across a range of values
#' This is useful to work wiht predict or add_epred_draw in a bayesian context
#' 
#' If you are want the marginal effects, consideration as to how to align
#' different predictor vectors should be considered
#' 
#' Already works well for categorical 
#' 
#' 
make_prediction_data <- function(data, 
                                 col_names = NULL,
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
    # hold_list <- list() # init list
    # 
    # for(col in names(char_cols)) {
    #   hold_list[[col]] <- unique(data[[col]])
    # }
    # #sequence is important here
    udf <- unique(data[,char_cols]) # make unique parings df
    # for len.out
    outdf <- udf[rep(1:nrow(udf), each = len.out),]

    # simple approach 
    if(!constrain) {
      if(nrow(udf) == nrow(data)) {
        warning('Character groups are all unique observations.')
      }
      num_df <- numdf_expand(data[,-char_cols], len.out = len.out, na.rm = T)
      num_df <- num_df[rep(1:nrow(num_df),times = nrow(udf)), ]
      
      outdf <- cbind(outdf, num_df)
      
    } else {
      
      if(nrow(udf) == nrow(data)) {
        stop('Character groups are all unique observations. Cannot constrain')
      }
      
      num_df <- as.data.frame(matrix(nrow = 0, ncol = ncol(data[,-char_cols])))
      names(num_df) <- names(data[,-char_cols])
      for(i in 1:nrow(udf)) {
        temp_df <- data #make temporary to filter from
        for(name in names(char_cols)) {
          temp_df <- temp_df |> 
            filter(.data[[name]] == udf[[name]][i])
        }
        num_df <- num_df |> 
          rbind(numdf_expand(temp_df[,-char_cols], len.out = len.out, na.rm = T))
      }
      
      outdf <- cbind(outdf,num_df)
      
      
    }
    
  } else {
    # this is simpler, there are no character
    outdf <- numdf_expand(data, len.out = len.out, na.rm = T)
  }
  
  return(outdf)
  
}

#' For all numeric df expansion

numdf_expand <- function(df, len.out = 1000, na.rm = T) {
  if(!all(sapply(df, is.numeric))) {
    stop('There are non-numeric data passed to numeric expansion')
  }
  var_names <- names(df)
  
  all_seqs <- var_names |> 
    lapply(function(name) {
      seq(from = min(df[[name]], na.rm = na.rm),
          to = max(df[[name]], na.rm = na.rm),
          length.out = len.out)
    }) |> 
    do.call(what = cbind,) |> 
    as.data.frame()
  
  names(all_seqs) = var_names
  return(all_seqs)
}


#' will not be perfect due to floating point issues i think
unscale <- function(vector, mu = NULL, sigma = NULL){
  
  if(is.null(mu)) {
    mu <- attr(vector, 'scaled:center')
    if(is.null(mu)){
      stop('No default mean found in provided vector')
    }
  }
  if(is.null(sigma)) {
    sigma <- attr(vector, 'scaled:scale')
    if(is.null(scale)) {
      stop('No default sd found in vector')
    }
  }
  
  y = (vector * sigma) + mu
  y = as.vector(y)
  return(y)
}
