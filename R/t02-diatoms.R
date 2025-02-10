rm(list = ls())
library(EcotaxaTools) # this is my package
library(dplyr)
library(ggplot2)
library(tidyr)
library(rstan)
library(lubridate)
source('./R/utils.R')

conc <- readRDS('./data/t01-copano_concentrations.rds')
name_list = readRDS('./data/t01-taxa_names.rds')
indv <- readRDS('./data/t01-copano_individals.rds')


##################
# MARK: Diatom Formating
##################

## \- Diatoms --------

all_diat <- conc |>
  filter(taxa %in% name_list$diatom) |> 
  group_by(sample_id) |> 
  summarize(count = sum(count), num_L = sum(num_L), um3_L = sum(um3_L)) |> 
  left_join(
    conc |> 
    select(-c(count, num_L, um3,um3_L, taxa)) |> 
      unique(),
    by = "sample_id"
  )

# impute to mean
all_diat$sal[which(is.na(all_diat$sal))] <- mean(all_diat$sal, na.rm = T)

# region \- simple plots

plot(all_diat$um3_L ~ all_diat$sal)
abline(lm(all_diat$um3_L ~ all_diat$sal))

plot(log10(all_diat$um3_L+1) ~ all_diat$sal)
abline(lm(log10(all_diat$um3_L+1) ~ all_diat$sal))

##################
# MARK: LogN Regression
##################



ln_mod_data <- list(
  N = nrow(all_diat[-which(is.na(all_diat$sal)),]),
  y = log10(all_diat$um3_L[-which(is.na(all_diat$sal))] + 1),
  sal = all_diat$sal[-which(is.na(all_diat$sal))]
)

all_diat$sal[is.na(all_diat$sal)] <- mean(all_diat$sal, na.rm = T)

ln_mod_data <- list(
  N = nrow(all_diat),
  y = log10(all_diat$um3_L + 1),
  sal = all_diat$sal
)

ln_fit <- stan(
  file = './stan/01-simple-sal_only.stan',
  data = ln_mod_data,
  iter = 2000,
  warmup = 500,
  cores = 14
)

# region \- pred & plot -----------

sal_pred = seq(min(all_diat$sal, na.rm = T), max(all_diat$sal, na.rm = T), by = 0.1)

ln_post <- post_predict(extract(ln_fit), list(s = sal_pred), "10^(beta0 + beta1 * s)") |> 
  summarize_pred()

ggplot() +
  geom_point(
    aes(
      x = all_diat$sal,
      y = log10(all_diat$um3_L + 1)
    )
  ) +
  geom_line(
    aes(
      x = sal_pred,
      y = log10(ln_post$mean)
    )
  ) +
  geom_ribbon(
    aes(
      x = sal_pred,
      ymin = log10(ln_post$low),
      ymax = log10(ln_post$high)
    )
  )

ln_rmse <- extract(ln_fit, pars = 'rmse_pred')
ln_rmse |> HDInterval::hdi()

###################
# MARK: Hurdle Model
##################

all_diat$z <- ifelse(all_diat$count > 0, 1, 0)


# region \- using BRMS ----

diat <- brm(
  bf(um3_L ~ sal,
     hu ~ sal),
  prior = c(
    set_prior('normal(0,1)', class = 'b', coef = 'sal', dpar = 'mu'),
    set_prior('normal(0,1)', class = 'b', dpar = 'mu'),
    set_prior('normal(0,1)', class = 'Intercept')
  ),
  family = hurdle_lognormal(),
  data = all_diat,
  thin = 2,
  chains = 3,
  iter = 999
)

# region \- Stan fit -------
hurdl_data <- list(
  N_obs = nrow(all_diat),
  z = ifelse(all_diat$count > 0, 1, 0),
  y = log10(all_diat$um3_L),
  sal = all_diat$sal  
)

hurd_fit <- stan(
  file = './stan/02-hurdle-sal_only.stan',
  data = hurdl_data,
  chains = 4, iter = 2000, warmup = 500, cores = 14
)


hurd_post <- post_predict(
  extract(hurd_fit),
  list(s = sal_pred),
  '(1/(1+ exp(alpha0 + alpha1 * s))) * 10^(beta0 + beta1 * s)'
) |> 
  summarize_pred()


ggplot() +
  geom_point(
    aes(
      x = all_diat$sal,
      y = all_diat$um3_L
    )
  ) +
  geom_line(
    aes(
      x = sal_pred,
      y = hurd_post$mean
    )
  ) +
  geom_ribbon(
    aes(
      x = sal_pred,
      ymin = hurd_post$low,
      ymax = hurd_post$high
    ),
    alpha = 0.5
  )


##################
# MARK: Mixture Approach
##################

indv_metrics <- indv |> 
  filter(taxo_name %in% name_list$diatom)

# region \- stan fit ----------
indv_metrics$sal[is.na(indv_metrics$sal)] <- mean(indv_metrics$sal, na.rm = T)
pois_data <- list(
  N_obs = nrow(all_diat),
  N_mes = nrow(indv_metrics),
  n = all_diat$count,
  img_vol = all_diat$img_vol/1000, # convert to L,
  sal = all_diat$sal,
  log_b = log(indv_metrics$um3),
  sal_b = indv_metrics$sal,
  temp_b = indv_metrics$t,
  obs_biomass_conc = all_diat$um3_L
)

pois_fit <- stan(
  file = './stan/03-mix-sal_only.stan',
  data = pois_data,
  chains = 4, iter = 5000, warmup = 500, cores = 14
)

# region \- pois plot

y_pred_c <- post_predict(
  extract(pois_fit),
  list(
    sal = sal_pred
  ),
  'exp(beta0 + beta1 * sal)'
) |> 
  summarize_pred()


ggplot() +
  geom_point(
    data = all_diat[all_diat$num_L < max(all_diat$num_L),],
    aes(
      x = sal,
      y = num_L
    ),
    size = 2.5,
    alpha = 0.5
  ) +
  geom_line(
    aes(
      x = sal_pred,
      y = y_pred_p$mean
    ),
    color = 'green',
    linewidth = 2
  )+
    geom_ribbon(
      aes(
        x = sal_pred,
        ymin = y_pred_c$low,
        ymax = y_pred_c$high
      ),
      fill = 'green',
      alpha = 0.25
    ) +
    labs(x = "Salinity", y = "Num / L")+
    theme_bw()

y_pred_b <- post_predict(
  extract(pois_fit),
  list(
    sal = sal_pred
  ),
  'exp(beta0 + beta1 * sal) * exp(a0)'
) |> 
  summarize_pred()


# region \- complete plot
ggplot() +
  geom_point(
    data = ,
    aes(
      x = all_diat$sal,
      y = all_diat$um3_L
    ),
    size = 2.5,
    alpha = 0.5
  )+
  geom_line(
    aes(
      x = sal_pred,
      y = y_pred_b$mean
    ),
    color = 'green',
    linewidth = 1
  )+
    geom_ribbon(
      aes(
        x = sal_pred,
        ymin = y_pred_b$low,
        ymax = y_pred_b$high
      ),
      fill = 'green',
      alpha = 0.25
    ) +
    # geom_line(
    #   aes(
    #     x = sal_pred,
    #     y = hurd_post$mean
    #   )
    # ) +
    # geom_ribbon(
    #   aes(
    #     x = sal_pred,
    #     ymin = hurd_post$low,
    #     ymax = hurd_post$high
    #   ),
    #   alpha = 0.5
    # )+
    labs(x = "Salinity", y = "um3 / L")+
    theme_bw()


ggplot() +
  geom_point(
    data = ,
    aes(
      x = all_diat$sal,
      y = log(all_diat$um3_L+1)
    ),
    size = 2.5,
    alpha = 0.5
  )+
  geom_line(
    aes(
      x = sal_pred,
      y = log(y_pred_b$mean)
    ),
    color = 'green',
    linewidth = 1
  )+
    geom_ribbon(
      aes(
        x = sal_pred,
        ymin = log(y_pred_b$low),
        ymax = log(y_pred_b$high)
      ),
      fill = 'green',
      alpha = 0.25
    ) +
    geom_line(
      aes(
        x = sal_pred,
        y = log(hurd_post$mean)
      )
    ) +
    geom_ribbon(
      aes(
        x = sal_pred,
        ymin = log(hurd_post$low),
        ymax = log(hurd_post$high)
      ),
      alpha = 0.5
    )+
    labs(x = "Salinity", y = "um3 / L")+
    theme_bw()
 