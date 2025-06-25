rm(list= ls())
library(cmdstanr)
library(dplyr)
library(tidyr)
library(ggplot2)

source('./R/utils-mod_making.R')


etx <- readRDS("./data/02-full_merged.RDS")


###################
# MARK: MODEL CONSTRUCTION 
###################



# region \-\- all together -------------
 
main_model <- cmdstan_model('./stan/gen-hierachical_mixture_model.stan')

bmass_priors <- list(
  b_mu = 8,
  b_sig = 3,
  tau_mu = 1,
  tau_sig = 1.5
)


# region \- mz model ------------
mz_poss_preds <- c('sal','t','do','turb','CHLA_N','windspeed', "TotalPAR")

mz_data <- prep_component_data('mz', mz_poss_preds)
mz_mod <- fit_mod(
  data_list = mz_data,
  count_preds = mz_poss_preds,
  bmass_priors,
  main_model,
  iter = 4000,
  chains = 5,
  parallel_chains = 5,
  iter_warmup = 2000,
  thin = 5,
  refresh = 500
)

# check pval 
pdf <- mz_mod$draws(c('MSE_obs_n','MSE_mod_n'), format = 'df')
p <- mean(pdf$MSE_obs_n > pdf$MSE_mod_n)
pdfb <- mz_mod$draws(c('MSE_obs_b','MSE_mod_b'), format = 'df')
pb <- mean(pdfb$MSE_obs_b > pdfb$MSE_mod_b)


# region \-\- check effect -------------
# note this checks the currrent model not candidate!
betas <- mz_mod$draws('beta_count', format = 'df')
theta_b <- mz_mod$draws('theta_bio', format = 'df')

saveRDS(
  list(
    beta = betas,
    mu_b = theta_b
  ),
  './data/03-post_mz.rds'
)
# region 