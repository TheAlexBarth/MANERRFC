rm(list= ls())
library(rstan)
library(dplyr)
library(tidyr)
library(ggplot2)

source('./R/utils-mod_making.R')


etx <- readRDS("./data/02-full_merged.RDS")

taxa_group <- 'diatom'
niter = 10000
nburn = 1000


core_data <- prep_component_data(taxa_group)

###################
# MARK: MODEL CONSTRUCTION 
###################

# region \- DEFINE predictors ------------
counts_poss_preds <-
indv_poss_preds <- 


# region \- possible preds --------------------------------

# region \-\- all together -------------
all_preds <- c(counts_poss_preds, counts_poss_preds)

 
main_model <- stan_model(file = './stan/gen-hierachical_mixture_model.stan')

mod_fit <- fit_component_mod(
  all_preds,
  core_data,
  main_model,
  chains = 1,
  iter = niter,
  warmup = nburn,
  cores = parallel::detectCores()
)

saveRDS(mod_fit, paste0('04-gen_mod_fit-', taxa_group,".RDS"))

########
# MARK: Diagnositcs
########

