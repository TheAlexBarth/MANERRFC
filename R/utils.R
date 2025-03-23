# MARK: Plot/Bonus

# region \- colorbline ggpalette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}



#' Colorblind friendly pallete
#' 
#' Define a color-blind friendly palette for up to 8 items.
#' Colors are from the
#' \href{http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette}{r-cookbook}
#' More than 8 will be defined by ggplot defaults.
#' 
#' @param n number to ID
#' 
#' @examples 
#' gg_cbb_col(3)
#' 
#' @export
gg_cbb_col <- function(n) {
  
  cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                  "#0072B2", "#D55E00", "#CC79A7")
  
  if(n > 8) {
    cbbPalette <- c(cbbPalette,gg_color_hue(n-8))
  }
  return(cbbPalette[1:n])
}

# MARK: Posterior Tools

# region \- prediction data 

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

#' region \- numeric expand

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

# MARK: Prediction from posterior 
post_predict <- function(posterior_draws, predict_data, formula, n_draws = NULL) {
  
  if(is.character(formula)) {
    formula <- parse(text = formula)
  }
  if(!inherits(formula, 'expression')) {
    stop('Formula must be an expression or string')
  }
  
  
  # draw length
  if(is.null(n_draws)) {
    n_draws <- posterior_draws |> 
      sapply(function(x) length(x)) |> 
      min()
  }
  
  
  # predictor length
  pred_length <- sapply(predict_data, length)
  if(length(unique(pred_length)) > 1) {
    stop('Predictor data must be the same length')
  }
  n_pred <- pred_length[1]
  
  # create environment for evaluation
  env <- c(posterior_draws, predict_data)
  
  pred <- matrix(NA, nrow = n_draws, ncol = n_pred)
  pred <- sapply(seq_len(n_pred), function(i) {
    env_updated <- env
    for(var in names(predict_data)) {
      env_updated[[var]] <- predict_data[[var]][i]
    }
    
    with(env_updated, eval(formula))
  })
  
  
  return(pred)
}

# region \- predition summary
# input a matrix of predition with columns as values, rows as draws
summarize_pred <- function(mat, detailed = TRUE, quantile = TRUE) {
  require(bayestestR)
  
  if(!detailed) {
    bounds <- apply(mat, 2, bayestestR::hdi)
    odf <- data.frame(
      mean = apply(mat, 2, mean),
      low = sapply(bounds, function(x) x$CI_Low[1]),
      high = sapply(bounds, function(x) x$CI_High[1])
    )
  } else if(quantile){
    qs = apply(mat, 2, quantile,probs = c(0.025,0.125,0.25,0.75,0.875,0.975))
    odf <- data.frame(
      median = apply(mat, 2, median),
      mean = apply(mat, 2, mean),
      low.50 = qs[3,],
      high.50 = qs[4,],
      low.75 = qs[2,],
      high.75 = qs[5,],
      low.95 = qs[1,],
      high.95 = qs[6,]
    )
  } else {
    bounds <- apply(mat, 2, bayestestR::hdi, ci = c(0.51,.75,.95))
    odf <- data.frame(
      median = apply(mat, 2, median),
      mean = apply(mat, 2, mean),
      low.50 = sapply(bounds, function(x) x$CI_low[1]),
      high.50 = sapply(bounds, function(x) x$CI_high[1]),
      low.75 = sapply(bounds, function(x) x$CI_low[2]),
      high.75 = sapply(bounds, function(x) x$CI_high[2]),
      low.95 = sapply(bounds, function(x) x$CI_low[3]),
      high.95 = sapply(bounds, function(x) x$CI_high[3])
    )
  }
  return(odf)
  
}

# region \- plot prediction ------------------

geom_error_range <- function(x, df, color) {
  out <- list(
    geom_ribbon(
      aes(
        x = x,
        ymin = df$low.95,
        ymax = df$high.95
      ),
      fill = color, alpha = 0.25
    ),
    geom_ribbon(
      aes(
        x = x,
        ymin = df$low.75,
        ymax = df$high.75
      ),
      fill = color, alpha = 0.25
    ),
    geom_ribbon(
      aes(
        x = x,
        ymin = df$low.50,
        ymax = df$high.50
      ),
      fill = color, alpha = 0.25
    )
  )
  return(out)
}

# region \- prediction list ---------
sum_list <- function(mat_list, sp_range = sal_pred){
  l = mat_list |> 
    lapply(summarize_pred) |> 
    lapply(
      function(x) mutate(x, sal = sp_range)
    ) |> 
    list_to_tib()|> 
    group_by(sal) |> 
    summarize(
      mean = sum(mean),
      low.50 = sum(low.50),
      high.50 = sum(high.50),
      low.75 = sum(low.75),
      high.75 = sum(high.75),
      low.95 = sum(low.95),
      high.95 = sum(high.95)
    ) |> 
    ungroup()
  
  return(l)
}

mean_list <- function(mat_list, sp_range = sal_pred){
  l = mat_list |> 
    lapply(summarize_pred) |> 
    lapply(
      function(x) mutate(x, sal = sp_range)
    ) |> 
    list_to_tib()|> 
    group_by(sal) |> 
    summarize(
      mean = mean(mean),
      low.50 = mean(low.50),
      high.50 = mean(high.50),
      low.75 = mean(low.75),
      high.75 = mean(high.75),
      low.95 = mean(low.95),
      high.95 = mean(high.95)
    ) |> 
    ungroup()
  
  return(l)
}