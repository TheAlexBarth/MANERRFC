#########
# Utilities to fit component models
#########

############
# MARK: Data prep tools 
###########



# region \- constr dfs --------------------------------
#' Function to prepare data for components
#' 
#' Needs taxonomic group and etx to exist in global scope
prep_component_data <- function(taxo_group) {
  
  # some warning
  if(!exists("etx")) {
    stop('etx does not exist or is incorrectly named')
  }
  if(!(taxo_group %in% names(etx$names))) {
    stop('taxo_groupo variable not correctly named with etx')
  }

  # region \-\- prep concentration --------------------
  conc_df <- etx$conc |> 
    filter(taxa %in% etx$names[[taxo_group]])
 
  # flip taxa names 
  conc_df$group <- conc_df$taxa |> 
    sapply(
      function(x) names(etx$names[[taxo_group]])[which(etx$names[[taxo_group]] == x)]
    ) |> 
    as.factor()

  all_conc <- conc_df |>
    group_by(sample_id, group, sample_site) |> 
    reframe(
      count = sum(count), 
      num_L = sum(num_L), 
      pgC_L = sum(pgC_L),
      img_vol = unique(img_vol),
      sal = unique(sal),
      t = unique(t),
      do = unique(DO),
      turb = unique(Turb),
      P = unique(P),
      NH4 = unique(NH4),
      N = unique(N),
      windspeed = unique(windspeed)
    )
  
  par_names <- c('sal', 't', 'do','turb','P','NH4','N','windspeed')

  # impute to mean
  for(par in par_names) {
    if(any(is.na(all_conc[[par]]))) {
      all_conc[[par]][which(is.na(all_conc[[par]]))] <- mean(all_conc[[par]], na.rm = T)
    }
  }

  # region \-\- indvidiual df -----------------------
  indv_mes <- etx$indv |> 
  filter(taxo_name %in% etx$names[[taxo_group]])

  indv_mes$group <- indv_mes$taxo_name |> 
    sapply(
      function(x) names(etx$names[[taxo_group]])[which(etx$names[[taxo_group]] == x)]
    ) |> 
    as.factor()

  for(par in par_names) {
    if(par %in% names(indv_mes)) {
      if(any(is.na(indv_mes[[par]]))) {
        indv_mes[[par]][which(is.na(indv_mes[[par]]))] <- mean(indv_mes[[par]], na.rm = T)
      }
    }
  }

  if(any(levels(all_conc$group) != levels(indv_mes$group))) {
    stop('There is something wrong with the levels!')
  } else {
    return(
      list(
        all_conc = all_conc,
        indv_mes = indv_mes
      )
    )
  }
}


# region \- make predictor matrices -------
#' Function to create predictor matrices
#' 
#' requires vectors for predictor names
#' AND prep_component_data list
prep_pred_matrix <- function(count_preds, indv_pred, data_list) {
  x_c <- 1
  for(pred in count_preds) {
    x_c <- cbind(x_c, data_list$all_conc[[pred]] |> scale())
  }
  colnames(x_c) = c('int',count_preds)

  x_i <- 1
  for(pred in indv_pred) {
    x_i <- cbind(x_i, data_list$indv_mes[[pred]] |> scale())
  }
  colnames(x_i) = c('int',indv_pred)

  pred_indv <- data_list$all_conc |> 
    ungroup() |> 
    select(all_of(indv_pred)) |> 
    mutate(
      across(everything(), scale)
    )
  if(any(is.na(pred_indv))){
    is.na(pred_indv) <- 0 # since it's scaled mean should be near 0
  }
  pred_indv <- cbind(1, pred_indv)

  return(
    list(
      X_counts = x_c,
      X_indv = x_i,
      pred_indv = pred_indv
    )
  )
}


###########
# MARK: model Fitting
###########

#' main model fitting loop
#' 
#' REQUIRES:
#' Pred list from prep_pred_matrix
#' data list from prep_component_data
#' compiled model
#' ... for chains, iters, cores, warmup
fit_component_mod <- function(pred_list, data_list, cmpld_model, ...) {

  stan_data <- list(
    N_obs = nrow(data_list$all_conc),
    N_mes = nrow(data_list$indv_mes),
    N_groups = length(levels(data_list$all_conc$group)),
    n = as.integer(data_list$all_conc$count),
    log_b = log(data_list$indv_mes$cmass),
    K_count = ncol(pred_list$X_counts),
    K_bmass = ncol(pred_list$X_indv),
    img_vol = data_list$all_conc$img_vol/1000,
    X_count = pred_list$X_counts,
    X_bio = pred_list$X_indv,
    group_counts = as.numeric(data_list$all_conc$group),
    group_bio = as.numeric(data_list$indv_mes$group),
    y_real = data_list$all_conc$pgC_L,
    Xpred_indv = pred_list$pred_indv
  )
  
  mod_fit <- sampling(
    object = cmpld_model,
    data = stan_data,
    ...
  )
  
  return(mod_fit)
}
