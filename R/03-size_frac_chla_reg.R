# regression for chla
#
# Predictors are raw nutrient concentrations (P, NH4, N, SiOH), not molar
# ratios - nutrient stoichiometry (N:P, N:Si) is reported separately in
# fig-02-environ_ts.R rather than used as a regression covariate.

rm(list = ls())
library(cmdstanr)
library(posterior)
library(dplyr)
library(lubridate)
library(tidyr)
source('./R/utils.R')


wq <- readRDS('./data/01a-swmp_wq_data.rds')
chl <- readRDS('./data/01b-size_frac_chl.RDS')
wind <- readRDS('./data/01c-wind_score.RDS')
si <- readRDS('./data/01e-silicate.rds')


all_data <- chl |>
  select(-year,-date) |>
  left_join(
    wq,
    by = c('yearmo','site'='sampling_site')
  ) |>
  left_join(
    wind |>
      select(yearmo, wind = wind_pca),
    by = 'yearmo',
    relationship = 'many-to-many'
  ) |>
  left_join(
    si |>
      select(SiOH, sampling_site, yearmo),
    by = c('yearmo','site'='sampling_site')
  )

all_data$month <- month(all_data$yearmo)

# define seasonal terms
all_data$sin_term <- sin(2*pi*all_data$month/12)
all_data$cos_term <- cos(2*pi*all_data$month/12)

all_data <- all_data |>
  pivot_longer(
    cols = c(micro, nano, pico),
    names_to = 'frac',
    values_to = 'ugL_chl'
  ) |>
  filter(!is.na(ugL_chl))


resp <- all_data$ugL_chl
frac <- all_data$frac |> factor(levels = levels(size_factors))
site <- all_data$site |> factor(levels = levels(site_factors))
preds <- all_data |>
  select(sin_term, cos_term, wind,temp, sal, P, NH4, N, SiOH)

# # log-transform nutrient terms
preds$P <- log(preds$P + 1e-5)
preds$NH4 <- log(preds$NH4 + 1e-5)
preds$N <- log(preds$N + 1e-5)
preds$SiOH <- log(preds$SiOH + 1e-5)

pred_scaled <- scale(preds[,-c(1,2,3)]) # don't scale seasonal or wind (already normalized)



# region Stan fit -----------------------
chl_model <- cmdstan_model('./stan/chl-norm_model.stan')

chl_data <- list(
  N_obs = length(resp),
  K_frac = length(levels(frac)),
  K_site = length(levels(site)),
  K_x = ncol(preds) + 1,
  y = resp,
  frac = as.numeric(frac),
  site = as.numeric(site),
  X = cbind(1, preds[,c(1,2,3)], pred_scaled) # pull back in seasonal terms
)

init_list <- list(
  sigma_f = rep(1, chl_data$K_frac),
  mu = matrix(0, nrow = chl_data$K_frac, ncol = chl_data$K_x),
  sigma_beta = matrix(1, nrow = chl_data$K_frac, ncol = chl_data$K_x),
  beta = array(0, dim = c(chl_data$K_frac, chl_data$K_site, chl_data$K_x))
)

niter= 5000
nthin = 5
nchain = 5
fit <- chl_model$sample(
  data = chl_data,
  init = rep(list(init_list),5),
  iter = niter,
  iter_warmup = 2500,
  chains = nchain,
  parallel_chains = nchain,
  thin = nthin,
  refresh = 100
)

# endregion -----------------------

# region posterior format -----------------------

# model check
dev <- fit$draws(c("dev_obs",'dev_sim'), format = 'df')
p_val <- dev$dev_obs > dev$dev_sim
cat('\n', paste0('Bayes pvalue: ', mean(p_val)), '\n')

cat('\n', 'Rhat summary:','\n',fit$summary()$rhat |> summary(), '\n')

betas <- fit$draws('beta')

beta_arr <- array(as.vector(betas), dim = c(niter/nthin * nchain, chl_data$K_frac, chl_data$K_site, chl_data$K_x))

sigma_f <- fit$draws('sigma_f')
sigma_arr <- array(as.vector(sigma_f), dim = c(niter/nthin * nchain, chl_data$K_frac))

saveRDS(
  list(
    beta = beta_arr,
    sigma_f = sigma_arr,
    data = list(
      all_data = all_data,
      X_scaled = cbind(1, preds[,c(1,2,3)], pred_scaled)
    )
  ),
  './data/03-post_chl.RDS'
)
