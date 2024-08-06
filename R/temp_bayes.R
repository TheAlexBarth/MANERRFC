####
# Playing around with brms hurdle modles #####
####


library(brms)
library(tidybayes)
library(emmeans)
library(broom)


###
# Chl ~ Sal #########
####


# basic linear regression with BRMS

# \- Model def ----------

all_data$chl_scaled <- scale(all_data$Chl)
all_data$sal_scaled <- scale(all_data$Sal)

chl_mod <- brm(
  chl_scaled ~ sal_scaled + sample_site + sample_site * sal_scaled,
  data = all_data,
  iter = 999,
  chains = 3,
  thin = 2
)

# predict_data <- data.frame(
#   chl_scaled = seq(min(all_data$chl_scaled, na.rm = T), 
#                    max(all_data$chl_scaled, na.rm = T),
#                    length.out = 100) |> 
#     rep(times = 2),
#   sal_scaled = seq(min(all_data$sal_scaled, na.rm = T),
#                    max(all_data$sal_scaled, na.rm = T),
#                    length.out = 100) |> 
#     rep(times = 2),
#   sample_site = rep(c('CE','CW'), each = 200)
# )

# \- Predictions & Plots -----------

chl_pred <- make_prediction_data(ungroup(all_data),
                                 col_names = c('sal_scaled',
                                               'chl_scaled',
                                               'sample_site'))



chl_pred <- add_epred_draws(chl_pred, chl_mod, ndraws = 500, dpar = T)

chl_pred$Sal <- chl_pred$sal_scaled |> 
  unscale(mu = mean(all_data$Sal,na.rm = T),
          sigma = sd(all_data$Sal, na.rm = T))

chl_pred$Chl <- chl_pred$chl_scaled |> 
  unscale(mu = mean(all_data$Chl, na.rm = T),
          sigma = sd(all_data$Chl, na.rm = T))

chl_pred$pred_chl <- chl_pred$.epred |> 
  unscale(mu = mean(all_data$Chl, na.rm = T),
          sigma = sd(all_data$Chl, na.rm = T))

## \-\- Overall ---------

ggplot()+
  stat_lineribbon(data = chl_pred,
                  aes(x = Sal, y = pred_chl, fill = sample_site),
                  alpha = 0.5, .width = 0.95) +
  geom_point(data = all_data,
             aes(x = Sal, y = Chl, color = sample_site)) + 
  labs(x = "Salinity [ppt]", y = 'Chlorophyll [ug/L]')+
  theme_bw()

## \-\- Marginal Effects --------------

site_pred <- make_prediction_data(ungroup(all_data), 'sample_site') ## ISSUE WITH CODE

site_pred <- data.frame(
  sample_site = c('CE','CW'),
  sal_scaled = 0
) |> 
  add_epred_draws(chl_mod, dpar = T, ndraws = 375)

site_pred$pred_chl <- unscale(site_pred$.epred,
                              mean(all_data$Chl, na.rm = T),
                              sd(all_data$Chl, na.rm = T))

ggplot(site_pred)+
  geom_density(aes(pred_chl, fill = sample_site)) +
  labs(x = 'Predicted Chl', y = 'Posterior Density', fill = "") +
  theme_bw()


posterior_chl <- chl_mod |> 
  as_draws_df(ndraw = 375)

b_SalCE <- (posterior_chl$b_sal_scaled)

b_SalCW <- (posterior_chl$b_sal_scaled + posterior_chl$`b_sal_scaled:sample_siteCW`)

chl_margins <- data.frame(
  b = b_SalCE,
  site = 'CE'
) |> 
  rbind(
    data.frame(
      b = b_SalCW,
      site = 'CW'
    )
  )

ggplot(chl_margins) +
  geom_density(aes(b, fill = site)) + 
  labs(x = 'Slope of Salinity', y= 'Posterior Density', fill = '')+
  theme_bw()

##
# Taxa-Specific Hurdles #######
##


## \- Diatoms --------

just_diat <- all_data |> 
  ungroup() |> 
  filter(taxa == 'Bacillariophyta')

diat <- brm(
  bf(num_L ~ sal_scaled,
     hu ~ sal_scaled + scale(img_vol)),
  prior = c(
    set_prior('normal(0,1)', class = 'b', coef = 'sal_scaled', dpar = 'mu'),
    set_prior('normal(0,1)', class = 'b', dpar = 'mu'),
    set_prior('normal(0,1)', class = 'Intercept')
  ),
  family = hurdle_lognormal(),
  data = just_diat,
  thin = 2,
  chains = 3,
  iter = 999
)

## \-\- Full Case -------

diat_pred <- make_prediction_data(just_diat,
                                  'sal_scaled')

## 
diat_pred$img_vol <- mean(scale(just_diat$img_vol), na.rm = T)

diat_pred <- diat_pred |> 
  add_epred_draws(diat, dpar = T, ndraw = 375)

diat_pred$Sal <- diat_pred$sal_scaled |> 
  unscale(mean(all_data$Sal, na.rm = T),
          sd(all_data$Sal, na.rm = T))

ggplot() +
  stat_lineribbon(data = diat_pred,
                  aes(x = Sal, y = log(.epred)),
                  .width = 0.95, alpha = 0.4,
                  color = 'purple', fill = 'purple', alpha = 0.5) +
  geom_smooth(data = just_diat,
              aes(x = Sal, y = log(num_L + 1)),
              method = 'lm', color = 'black', fill = "grey",
              alpha = 0.75) + 
  geom_point(data = just_diat,
             aes(x = Sal, y = log(num_L + 1))) +
  labs(x = 'Salinity', y = 'Ln(Density) #/L', fill = "") +
  theme_bw()

## \-\- Just mu model ----------------


ggplot() +
  stat_lineribbon(data = diat_pred,
                  aes(x = Sal, y = log(exp(mu)+1)),
                  .width = 0.95, alpha = 0.4) +
  geom_point(data = just_diat,
             aes(x = Sal, y = log(num_L))) +
  labs(x = 'Salinity', y = '#/L', fill = "") +
  theme_bw()


## \-\- Detection model -----------------


## \-\-\- Salinity ---------------
diat_detect_pred <- make_prediction_data(just_diat,
                                         "sal_scaled")


diat_detect_pred <- diat_detect_pred[rep(1:nrow(diat_detect_pred)),] |> 
  as.data.frame()

diat_detect_pred$img_vol <- mean(scale(just_diat$img_vol))

names(diat_detect_pred) <- c('sal_scaled', 'img_vol')

diat_detect_pred <- diat_detect_pred |> 
  add_epred_draws(diat, ndraw = 300, dpar = T)

diat_detect_pred$Sal <- unscale(diat_detect_pred$sal_scaled,
                                mean(all_data$Sal, na.rm = T),
                                sd(all_data$Sal, na.rm = T))

diat_detect_pred$img_vol <- unscale(diat_detect_pred$img_vol,
                                    mean(just_diat$img_vol, na.rm = T),
                                    sd(just_diat$img_vol, na.rm = T))

ggplot() +
  stat_lineribbon(data = diat_detect_pred,
                  aes(x = Sal, y = 1-hu),
                  .width = 0.95, alpha = 0.4) +
  geom_point(data = just_diat,
             aes(x = Sal, y = as.numeric((num_L>0)))) +
  labs(x = 'Salinity', y = '#/L', fill = "") +
  theme_bw()


## \-\-\- Img Vol ---------------------------

just_diat$img_vol_s <- scale(just_diat$img_vol)
img_vol_predict <- make_prediction_data(just_diat,
                                        'img_vol_s')
  
img_vol_predict$sal_scaled <- mean(just_diat$sal_scaled, na.rm = T)

names(img_vol_predict) <- c('sal_scaled', 'img_vol') |> 
  rev()

img_vol_predict <- img_vol_predict |> 
  add_epred_draws(diat, ndraw = 300, dpar = T)


img_vol_predict$img_vol <- unscale(img_vol_predict$img_vol,
                                    mean(just_diat$img_vol, na.rm = T),
                                    sd(just_diat$img_vol, na.rm = T))


ggplot() +
  stat_lineribbon(data = img_vol_predict,
                  aes(x = img_vol, y = 1-hu),
                  .width = 0.95, alpha = 0.4) +
  geom_point(data = just_diat,
             aes(x = img_vol, y = as.numeric((num_L>0)))) +
  labs(x = 'Sampling Volum', y = '1-p', fill = "") +
  theme_bw()



# \- Dinos ------------------


just_dinos <- all_data |> 
  ungroup() |> 
  filter(taxa == 'Dinophyceae')

dinos <- brm(
  bf(num_L ~ sal_scaled,
     hu ~ sal_scaled + scale(img_vol)),
  family = hurdle_lognormal(),
  data = just_dinos,
  thin = 2,
  chains = 3,
  iter = 999
)

## \-\- Full Case -------

dinos_pred <- make_prediction_data(just_dinos,
                                  'sal_scaled')

## 
dinos_pred$img_vol <- mean(scale(just_dinos$img_vol), na.rm = T)

dinos_pred <- dinos_pred |> 
  add_epred_draws(dinos, dpar = T, ndraw = 375)

dinos_pred$Sal <- dinos_pred$sal_scaled |> 
  unscale(mean(all_data$Sal, na.rm = T),
          sd(all_data$Sal, na.rm = T))

ggplot() +
  stat_lineribbon(data = dinos_pred,
                  aes(x = Sal, y = log(.epred)),
                  .width = 0.95, alpha = 0.4) +
  stat_smooth(data = just_dinos,
              aes(x = Sal, y = log(num_L)),
              method = 'lm')+
  geom_point(data = just_dinos,
             aes(x = Sal, y = log(num_L))) +
  labs(x = 'Salinity', y = '#/L', fill = "") +
  theme_bw()

## \-\- Just mu model ----------------


ggplot() +
  stat_lineribbon(data = dinos_pred,
                  aes(x = Sal, y = mu),
                  .width = 0.95, alpha = 0.4) +
  geom_point(data = just_dinos,
             aes(x = Sal, y = log(num_L))) +
  labs(x = 'Salinity', y = '#/L', fill = "") +
  theme_bw()
