####
# Playing around with brms hurdle modles #####
####


library(brms)
library(tidybayes)
library(emmeans)
library(broom)


# basic linear regression with BRMS


all_data$chl_scaled <- scale(all_data$Chl)
all_data$sal_scaled <- scale(all_data$Sal)

chl_mod <- brm(
  chl_scaled ~ sal_scaled + sample_site + sal_scaled*sample_site,
  data = all_data,
  iter = 999,
  chains = 3,
  thin = 2
)

predict_data <- data.frame(
  chl_scaled = seq(min(all_data$chl_scaled, na.rm = T), 
                   max(all_data$chl_scaled, na.rm = T),
                   length.out = 100) |> 
    rep(times = 2),
  sal_scaled = seq(min(all_data$sal_scaled, na.rm = T),
                   max(all_data$sal_scaled, na.rm = T),
                   length.out = 100) |> 
    rep(times = 2),
  sample_site = rep(c('CE','CW'), each = 200)
)

test_pred <- epred_draws(chl_mod, predict_data, ndraws = 400)

ggplot(test_pred,
       aes(x = sal_scaled, y = .epred, fill = sample_site))+
  stat_lineribbon(alpha = 0.5, .width = 0.95)




model_boring <- brm(
  bf(lifeExp ~ year),
  data = gapminder,
  chains = 4,
  cores = 4
)

# 
# just_diat <- all_data |> 
#   filter(taxa == 'Bacillariophyta')
# 
# 
# just_dinos <- 
# 
# test <- brm(
#   bf(num_L ~ Sal + sample_site,
#      hu ~ img_vol + Sal),
#   family = hurdle_lognormal(),
#   prior = ()
#   data = just_diat,
#   thin = 2
# )
# 
# 
# emmeans()
