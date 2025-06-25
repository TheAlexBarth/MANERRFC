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
prep_component_data <- function(taxo_group, poss_preds) {
  require(dplyr)
  
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
      date = unique(date),
      img_vol = unique(img_vol),
      sal = unique(sal),
      t = unique(t),
      do = unique(DO),
      turb = unique(Turb),
      P = unique(P),
      NH4 = unique(NH4),
      N = unique(N),
      windspeed = unique(windspeed),
      TotalPAR = unique(TotalPAR),
      CHLA_N = unique(CHLA_N)
    )
  
  # region \-\- drop NA pars ----------------------
  drop_rows <- apply(
    all_conc, 1, function(x) {
      any(sapply(poss_preds, function(y) is.na(x[[y]])))
    }
  )

  all_conc <- all_conc[which(!drop_rows),]
  

  # region \-\- indvidiual df -----------------------
  indv_mes <- etx$indv |> 
  filter(taxo_name %in% etx$names[[taxo_group]])

  indv_mes$group <- indv_mes$taxo_name |> 
    sapply(
      function(x) names(etx$names[[taxo_group]])[which(etx$names[[taxo_group]] == x)]
    ) |> 
    as.factor()

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



# MARK: Site Specific --------------------------------
#' Function to prepare data for components
#' 
#' Needs taxonomic group and etx to exist in global scope
prep_component_data_sites <- function(taxo_group, poss_preds) {
  require(dplyr)
  
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
  
  conc_df$sample_site <- as.factor(conc_df$sample_site)

  all_conc <- conc_df |>
    group_by(sample_id, group, sample_site) |> 
    reframe(
      count = sum(count), 
      num_L = sum(num_L), 
      pgC_L = sum(pgC_L),
      date = unique(date),
      img_vol = unique(img_vol),
      sal = unique(sal),
      t = unique(t),
      do = unique(DO),
      turb = unique(Turb),
      P = unique(P),
      NH4 = unique(NH4),
      N = unique(N),
      windspeed = unique(windspeed),
      TotalPAR = unique(TotalPAR),
      CHLA_N = unique(CHLA_N)
    )
  
  # region \-\- drop NA pars ----------------------
  drop_rows <- apply(
    all_conc, 1, function(x) {
      any(sapply(poss_preds, function(y) is.na(x[[y]])))
    }
  )

  all_conc <- all_conc[which(!drop_rows),]
  

  # region \-\- indvidiual df -----------------------
  indv_mes <- etx$indv |> 
  filter(taxo_name %in% etx$names[[taxo_group]])

  indv_mes$group <- indv_mes$taxo_name |> 
    sapply(
      function(x) names(etx$names[[taxo_group]])[which(etx$names[[taxo_group]] == x)]
    ) |> 
    as.factor()

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
fit_mod <- function(data_list, count_preds, bmass_priors, cmpld_model, ...) {
  require(cmdstanr)

  x_pred <- 1
  for(pred in count_preds) {
    x_pred <- cbind(x_pred, scale(data_list$all_conc[[pred]]))
  }

  init_list <- list(
    theta_count = rep(0, length(count_preds)),
    tau_count = rep(0.5, length(count_preds)),
    theta_bio = 0,
    tau_bio = 0.5,
    beta_count = matrix(
      0, nrow = length(levels(data_list$all_conc$group)), ncol = length(count_preds)
    ),
    mu_bio = rep(bmass_priors$b_mu, length(levels(data_list$all_conc$group))),
    sigma_bio = 1
  )

  stan_data <- list(
    N_obs = nrow(data_list$all_conc),
    N_mes = nrow(data_list$indv_mes),
    N_groups = length(levels(data_list$all_conc$group)),
    n = as.integer(data_list$all_conc$count),
    log_b = log(data_list$indv_mes$cmass),
    K_count = ncol(x_pred),
    img_vol = data_list$all_conc$img_vol/1000, # convert to L,
    X_count = x_pred,
    group_counts = as.numeric(data_list$all_conc$group),
    group_bio = as.numeric(data_list$indv_mes$group),
    theta_b_mu = bmass_priors$b_mu,
    theta_b_sig = bmass_priors$b_sig,
    tau_b_mu = bmass_priors$tau_mu,
    tau_b_sig = bmass_priors$tau_sig
  )
  
  mod_fit <- cmpld_model$sample(
    data = stan_data,
    ...
  )
  
  return(mod_fit)
}

overlaps_zero <- function(mat) {
  apply(mat, 2, function(c){
    min(c) <= 0 && max(c) >= 0
  })
}
