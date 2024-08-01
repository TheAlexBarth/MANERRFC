set.seed(10)

just_diat

p_off <- brm(
  bf(count ~ sal_scaled + offset(log(img_vol))),
  prior = c(
    set_prior('normal(0,100)', class = 'Intercept'),
    set_prior('normal(0,1)', class = 'b')
  ),
  data = just_diat,
  family = poisson,
  iter = 999,
  chains = 3,
  thin = 2,
  cores = 4
)


poff_pred <- make_prediction_data(just_diat,
                                  'sal_scaled')

poff_pred$img_vol <- mean(just_diat$img_vol)

poff_pred <- poff_pred |> 
  add_epred_draws(p_off, ndraw = 400, dpar = T)

poff_pred$Sal <- unscale(poff_pred$sal_scaled,
                         mean(just_diat$Sal, na.rm = T),
                         sd(just_diat$Sal, na.rm = T))


ggplot() +
  stat_lineribbon(data = poff_pred,
                  aes(x = Sal, y = .epred),
                  .width = 0.95, alpha = 0.5) +
  geom_point(data = just_diat,
             aes(x = Sal, y = count)) +
  theme_bw()




### \-\- Compare models ---------

diat_loo <- loo(diat, save_psis = T) 

diat_den = posterior_predict(diat)

diat_plot <- bayesplot::ppc_loo_pit_qq(
  yrep = diat_den,
  y = just_diat$num_L[!is.na(just_diat$Sal)],
  lw = weights(diat_loo$psis_object)
)

poff_loo <- loo(p_off, save_psis = T)
poff_count <- posterior_predict(p_off)
poff_plot <- bayesplot::ppc_loo_pit_qq(
  yrep = poff_count,
  y = just_diat$num_L[!is.na(just_diat$Sal)],
  lw = weights(poff_loo$psis_object)
)
poff_plot

## \-\- Dinos --------

dino_poff <- brm(
    bf(count ~ sal_scaled + offset(log(img_vol)) + s_temp),
    prior = c(
      set_prior('normal(0,100)', class = 'Intercept'),
      set_prior('normal(0,1)', class = 'b')
    ),
    data = just_dinos,
    family = poisson,
    iter = 999,
    chains = 3,
    thin = 2,
    cores = 4
  )


just_dinos$s_temp <- scale(just_dinos$Temp)
poff_pred <- make_prediction_data(just_dinos,
                                  c('s_temp',
                                    'sal_scaled'))

poff_pred$img_vol <- mean(just_dinos$img_vol) 
poff_pred$sal_scaled <- 0

poff_pred <- poff_pred |> 
  add_epred_draws(dino_poff, ndraw = 400, dpar = T)

poff_pred$Sal <- unscale(poff_pred$sal_scaled,
                         mean(just_diat$Sal, na.rm = T),
                         sd(just_diat$Sal, na.rm = T))


ggplot() +
  stat_lineribbon(data = poff_pred,
                  aes(x = s_temp, y = .epred),
                  .width = 0.95, alpha = 0.5) +
  geom_point(data = just_dinos,
             aes(x = s_temp, y = count)) +
  theme_bw()





