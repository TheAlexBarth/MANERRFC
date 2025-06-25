

# region \- marginal plots ---------------------
make_marg_data <- function(param, data, mod, post_coef, sim.res = 1000) {
  require(EcotaxaTools)
  warning("I made this assuming no effects on biomass")
  # region \- safety check --------------------
  if(!(param %in% mod$count_preds)) {
    stop(paste0(param, ' is not available. Only can use: ', paste(mod$count_preds, collapse = ', ')))
  }
  if(!(all(names(post_coef) %in% c('v','beta_bio','beta_count')))) {
    stop('Something got fucked up with the posterior names.')
  }

  mean_vol <- mean(data$all_conc$img_vol/1000)

  # regoin \- format betas -------------------
  beta_list = list()
  for(i in 1:length(mod$group_names)) {
    beta_list[[i]] <- post_coef$beta_count[,i,]
    colnames(beta_list[[i]]) <- c('int',mod$count_preds) # guess unnecessary but ok
  }
  names(beta_list) <- mod$group_names

  # region \- simulate new range ------------------

  par_mean = mean(data$all_conc[[param]], na.rm = T)
  par_sd = sd(data$all_conc[[param]], na.rm = T)
  par_range = range(scale(data$all_conc[[param]]), na.rm = T)
  
  X_new <- matrix(0, nrow = sim.res, ncol = length(mod$count_preds))
  colnames(X_new) <- mod$count_preds
  X_new[, which(colnames(X_new) == param)] <- seq(
    par_range[1], par_range[2], length.out = sim.res
  )
  X_new <- cbind(1, X_new) # add intercept


  # region \- simulate loop ---------------------------

  g_list <- list() #group list
  lambda <- list() # lambda list
  n_list <- list() # count list
  # loop through all levels of group
  cat("", "=======", 'Simulating Counts', "========", "", sep = '\n')
  pb = txtProgressBar(0, length(mod$group_names))
  setTxtProgressBar(pb, 0)
  for(taxo in mod$group_names) {
    #make lambda matrix (draws of posterior estimates for param)
    lambda[[taxo]] <- exp(beta_list[[taxo]] %*% t(X_new) + log(mean_vol))

    #from lambda mat - simulate counts
    n_list[[taxo]] <- matrix(
      rpois(length(lambda[[taxo]]), lambda = lambda[[taxo]]),
      nrow = nrow(lambda[[taxo]]),
      ncol = ncol(lambda[[taxo]])
    )
    setTxtProgressBar(pb, which(taxo == mod$group_names))
  }


  cat("","=====", "Simulating Biomass","=====","", sep = '\n')
  setTxtProgressBar(pb, 0)
  for(i in 1:length(mod$group_names)) {
    taxo = mod$group_names[i]

    log_d <- rnorm(
      n = sim.res * nrow(lambda[[taxo]]),
      mean = matrix(post_coef$beta_bio[,i,]) %*% matrix(1, ncol = sim.res),
      sd = post_coef$v
    ) |> 
      matrix(nrow = nrow(lambda[[taxo]]), ncol = sim.res)

    d <- exp(log_d)

    # predict biomass from d and n.
    g_list[[taxo]] = (n_list[[taxo]] * d) / mean_vol
    setTxtProgressBar(pb, i)
  }

  # could add update to count sumammary if desired here...

  cat("","======","Finishing Up, slowly....", "======", "", sep = "\n")

  out_list <- g_list |> 
    lapply(summarize_pred) |> 
    lapply(
      function(x) mutate(x, "{param}" := unscale(X_new[,param], par_mean, par_sd))
    )
  
  out_list$total <- Reduce('+', g_list) |> 
    summarize_pred() |> 
    mutate("{param}" := unscale(X_new[,param], par_mean, par_sd))
  
  # this is slow
  # out_list$total = g_list |> 
  #   list_to_tib() |>
  #   group_by(across(all_of(param))) |> 
  #   reframe(
  #     mean = sum(mean),
  #     low.50 = sum(low.50),
  #     high.50 = sum(high.50),
  #     low.75 = sum(low.75),
  #     high.75 = sum(high.75),
  #     low.95 = sum(low.95),
  #     high.95 = sum(high.95)
  #   ) |> 
  #   ungroup()

  return(out_list)
}


make_marg_plot <- function(taxa, sim_data, fill_col, ...) {
  require(ggplot2)
  require(dplyr)

  par_name = names(sim_data[[taxa]])[9]

  if(!exists("group_data")) {
    stop('was group data deleted? make sure ran sequentially')
  }

  if(!(taxa %in% names(sim_data))) {
    stop('invalid taxa name')
  } else if(taxa != 'total') {
    filtered_data <- group_data$all_conc |> 
      filter(group == taxa)
  } else {
    filtered_data <- group_data$all_conc |> 
      group_by(across(all_of(c('sample_id',par_name)))) |> 
      reframe(
        pgC_L = sum(pgC_L)
      )
  }

  filtered_data <- filtered_data[filtered_data[[par_name]] != mean(filtered_data[[par_name]], na.rm = T),]
  p = ggplot() +
    geom_point(
      aes(
        x = filtered_data[[par_name]],
        y = filtered_data$pgC_L/5
      ),
      ...
    ) +
    geom_error_range(
      x = sim_data[[taxa]][[par_name]],
      df = sim_data[[taxa]],
      color = fill_col
    )
  return(p)
}


######
# MARK: TS DATA MAKER
######



make_ts_data <- function(site, wq, data, mod, post_coef) {
  
  require(EcotaxaTools)
  warning("I made this assuming no effects on biomass")
  # region \- safety check --------------------
  if(!(all(names(post_coef) %in% c('v','beta_bio','beta_count')))) {
    stop('Something got fucked up with the posterior names.')
  }

  mean_vol <- mean(data$all_conc$img_vol/1000)

  # regoin \- format betas -------------------
  beta_list = list()
  for(i in 1:length(mod$group_names)) {
    beta_list[[i]] <- post_coef$beta_count[,i,]
    colnames(beta_list[[i]]) <- c('int',mod$count_preds) # guess unnecessary but ok
  }
  names(beta_list) <- mod$group_names

  # region \- scale_new_range ------------------

  wq <- wq |> 
    filter(sample_site == site)

  X_new <- wq[,mod$count_preds]
  for(name in names(X_new)) {
    X_new[[name]] <- scale(X_new[[name]], 
      center = mean(data$all_conc[[name]], na.rm = T),
      scale = sd(data$all_conc[[name]], na.rm = T)
    )
  }
  X_new <- as.matrix(cbind(1, X_new)) # add intercept


  # region \- simulate loop ---------------------------

  g_list <- list() #group list
  lambda <- list() # lambda list
  n_list <- list() # count list
  # loop through all levels of group
  cat("", "=======", 'Simulating Counts', "========", "", sep = '\n')
  pb = txtProgressBar(0, length(mod$group_names))
  setTxtProgressBar(pb, 0)
  for(taxo in mod$group_names) {
    #make lambda matrix (draws of posterior estimates for param)
    lambda[[taxo]] <- exp(beta_list[[taxo]] %*% t(X_new) + log(mean_vol))

    #from lambda mat - simulate counts
    n_list[[taxo]] <- matrix(
      rpois(length(lambda[[taxo]]), lambda = lambda[[taxo]]),
      nrow = nrow(lambda[[taxo]]),
      ncol = ncol(lambda[[taxo]])
    )
    setTxtProgressBar(pb, which(taxo == mod$group_names))
  }


  cat("","=====", "Simulating Biomass","=====","", sep = '\n')
  setTxtProgressBar(pb, 0)
  for(i in 1:length(mod$group_names)) {
    taxo = mod$group_names[i]

    log_d <- rnorm(
      n = nrow(wq) * nrow(lambda[[taxo]]),
      mean = matrix(post_coef$beta_bio[,i,]) %*% matrix(1, ncol = nrow(wq)),
      sd = post_coef$v
    ) |> 
      matrix(nrow = nrow(lambda[[taxo]]), ncol = nrow(wq))

    d <- exp(log_d)

    # predict biomass from d and n.
    g_list[[taxo]] = (n_list[[taxo]] * d) / mean_vol
    setTxtProgressBar(pb, i)
  }

  # could add update to count sumammary if desired here...

  cat("","======","Finishing Up, slowly....", "======", "", sep = "\n")

  out_list <- g_list |> 
    lapply(summarize_pred) |> 
    lapply(
      function(x) mutate(x, date = wq$date)
    )
  
  out_list$total <- out_list |> 
    do.call(what = rbind,) |> 
    group_by(date) |> 
    reframe(
      across(everything(), sum, na.rm = T)
    )
  
  # this is slow
  # out_list$total = g_list |> 
  #   list_to_tib() |>
  #   group_by(across(all_of(param))) |> 
  #   reframe(
  #     mean = sum(mean),
  #     low.50 = sum(low.50),
  #     high.50 = sum(high.50),
  #     low.75 = sum(low.75),
  #     high.75 = sum(high.75),
  #     low.95 = sum(low.95),
  #     high.95 = sum(high.95)
  #   ) |> 
  #   ungroup()

  return(out_list)
}



# need to specify taxa AND site
make_ts_plot <- function(taxa, site, sim_data, fill_col, ...) {
  require(ggplot2)
  require(dplyr)


  group_data$all_conc$date = group_data$all_conc$sample_id |> 
    sapply(function(x) gsub('[A-Z][A-Z]_', "",x)) |> 
    as.Date(format = '%Y-%m-%d')

  if(!exists("group_data")) {
    stop('was group data deleted? make sure ran sequentially')
  }

  if(!(taxa %in% names(sim_data))) {
    stop('invalid taxa name')
  } else if(taxa != 'total') {
    filtered_data <- group_data$all_conc |> 
      filter(group == taxa, sample_site == site)
  } else {
    filtered_data <- group_data$all_conc |> 
      filter(sample_site == site) |> 
      group_by(date) |> 
      summarize(
        pgC_L = sum(pgC_L)
      )
  }

  p = ggplot() +
    geom_point(
      aes(
        x = filtered_data$date,
        y = filtered_data$pgC_L
      ),
      ...
    ) +
    geom_error_range(
      x = sim_data[[taxa]]$date,
      df = sim_data[[taxa]],
      color = fill_col
    )
  return(p)
}
