rm(list = ls())
library(cmdstanr)
library(posterior)
library(dplyr)
library(lubridate)
library(tidyr)
source('./R/utils.R')
set.seed(20250630)

etx <- readRDS('./data/02-full_merged.rds')


# region Data prep -----------------------
conc <- etx$conc

conc$month <- conc$yearmo |> month()
conc$sin_term <- sin(2*pi*conc$month/12)
conc$cos_term <- sin(2*pi*conc$month/12)



preds <- conc |> 
  select(
    sin_term, cos_term, wind_pca,
    temp, sal, P, NH4, N,
    micro, nano, pico
  )

preds$P <- log(preds$P + 1e-5)
preds$NH4 <- log(preds$NH4 + 1e-5)
preds$N <- log(preds$N + 1e-5)

pred_scaled <- scale(preds[,-c(1,2,3)]) # don't scale seasonal or wind (already normalized)

# indexing group factors

role <- conc$functional_role
site <- conc$sample_site
taxo <- etx$indv$taxo_name |> 
  as.factor()

taxo_role <- data.frame(taxo_name = levels(taxo)) |> 
  left_join(
    etx$indv |> 
      select(taxo_name, functional_role) |> 
      distinct(),
    by = 'taxo_name',
    relationship = 'one-to-one'
  )

taxo_role <- factor(taxo_role$functional_role, levels = troph_factors)


incld_chl <- ifelse(grepl('_auto',conc$functional_role),0,1)

# endregion -----------------------


# region Stan fit -----------------------
two_mod <- cmdstan_model('./stan/site-two_part_mod.stan')


micro_data <- list(
  N_obs = length(conc$count),
  N_mes = length(etx$indv$cmass),
  K_role = length(levels(role)),
  K_taxo = length(levels(taxo)),
  K_site = length(levels(site)),
  K_x = ncol(preds) +1,
  #data
  n = conc$count |> as.integer(),
  img_vol = conc$img_vol, # ml
  W = cbind(1, preds[,c(1,2,3)], pred_scaled),
  b = etx$indv$cmass,
  # index
  role = role,
  site = site,
  incld_chl = incld_chl,
  taxo_role = taxo_role,
  taxo = taxo
)

init_list <- function() {
  list(
    gamma = rep(1, micro_data$K_role),
    alpha = array(
      rnorm(micro_data$K_role * micro_data$K_site* micro_data$K_x, 0, 0.1),
      dim = c(micro_data$K_role, micro_data$K_site, micro_data$K_x)
    ),
    mu_alpha = matrix(rnorm(micro_data$K_role * micro_data$K_x, 0,0.1), micro_data$K_role, micro_data$K_x),
    sigma_alpha = matrix(rexp(micro_data$K_role * micro_data$K_x), micro_data$K_role, micro_data$K_x),
    sigma_g = rexp(micro_data$K_taxo, 1),
    eta_g = rnorm(micro_data$K_taxo, 10^4, 10^2),
    eta_r = rnorm(micro_data$K_role, 10^4, 10^2),
    sigma_r = rexp(micro_data$K_role, 0.5)
  )  
}

niter = 5000
nthin = 5
nchain = 5

fit <- two_mod$sample(
  data = micro_data,
  init = lapply(1:nchain, function(x) init_list()),
  iter = niter,
  iter_warmup = niter/2,
  chains = nchain,
  parallel_chains = nchain,
  thin = nthin,
  refresh = 100
)

# endregion -----------------------

# region Posterior format -----------------------

mse <- fit$draws(c('MSE_obs','MSE_mod'), format = 'df')
pval_pois <- mse$MSE_obs > mse$MSE_mod
cat('\n', paste0('Bayes pvalue for regression: ', mean(pval_pois)), '\n') 

dev <- fit$draws(c('dev_obs','dev_mod'), format = 'df')
pval_bmass <- dev$dev_obs > dev$dev_mod
cat('\n', paste0('Bayes pvalue for biomass: ', mean(pval_bmass)), '\n')

cat('\n', 'Rhat summary:','\n',fit$summary()$rhat |> summary(), '\n') 



alphas <- fit$draws('alpha')
alpha_arr <- array(
  as.vector(alphas), 
  dim = c(niter/nthin * nchain, micro_data$K_role, micro_data$K_site, micro_data$K_x),
)

eta_r <- fit$draws('eta_r')
eta_mat <- array(as.vector(eta_r), dim = c(niter/nthin * nchain, micro_data$K_role))
sigma_r <- fit$draws('sigma_r')
sigma_mat <- array(as.vector(sigma_r), dim = c(niter/nthin*nchain, micro_data$K_role))

saveRDS(
  list(
    alpha = alpha_arr,
    eta_r = eta_mat,
    sigma_r = sigma_mat,
    preds = cbind(1, preds[,c(1,2,3)], pred_scaled)
  ),
  './data/04-post_micro.RDS'
)

# endregion -----------------------
