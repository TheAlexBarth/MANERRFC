rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

source("./R/utils.r")

# region Prep data -----------------------
post_chl <- readRDS('./data/03-post_chl.RDS')


# endregion -----------------------



# region Create effect table -----------------------

amp = sqrt(post_chl$beta[,,,2]^2 + post_chl$beta[,,,3]^2)
beta_main <- post_chl$beta[,,,-c(1,2)] # drop sin term
beta_main[,,,1] <- amp # replace cos with amp



group_levels = list(
  'frac' = levels(size_factors), 
  'site' = levels(site_factors), 
  'x' = c('seasonal', names(post_chl$data$X_scaled[,-c(1,2,3)]))
)

beta_df <- post_arr_to_df(beta_main, c('frac','site','x'), dim_levels = group_levels)


beta_sum <- beta_df |> 
  group_by(frac, site, x) |> 
  summarize(
    low = quantile(est, probs = c(0.025)),
    mid = median(est),
    high = quantile(est, probs = c(0.975)),
    pos_prob = mean(est > 0)
  )


effect_table <- beta_sum |> 
  select(x, mid, low, high, pos_prob) |> 
  mutate(
    `Post. Mean` = sprintf('%.2f', mid),
    `95% Cred.Intv` = paste0("(", sprintf("%.2f", low), ", ", sprintf("%.2f", high), ")"),
    `Pr(>0)` = sprintf("%.2f", pos_prob)
  ) |> 
  select(`Size Frac` = frac, `Site` = site, Variable = x, `Post. Mean`, `95% Cred.Intv`, `Pr(>0)`)

write.csv(effect_table, './output/s-chl_post_summary.csv', row.names = FALSE)

# endregion -----------------------

# region Time Series -----------------------

ts_list <- list()
for(frac in size_factors) {
  ts_list[[frac]] <- list()
  frac_idx <- as.numeric(size_factors[which(size_factors == frac)])
  for(site in site_factors) {
    site_idx <- as.numeric(site_factors[which(site_factors == site)])

    data_idx <- which(post_chl$data$all_data$frac == frac & post_chl$data$all_data$site == site)

    pred_chl <- as.matrix(
      post_chl$beta[,frac_idx, site_idx,]) %*% t(as.matrix(post_chl$data$X_scaled[data_idx,])
    )

    sim_chl <- matrix(
      rnorm(length(pred_chl), mean = pred_chl, sd = post_chl$sigma_f[,frac_idx]),
      ncol = ncol(pred_chl),
      nrow = nrow(pred_chl)
    ) |> 
      summarize_pred()
    sim_chl$yearmo = post_chl$data$all_data$yearmo[data_idx]
    ts_list[[frac]][[site]] <- sim_chl
  }
}


ggplot() +
  geom_error_range(
    x = ts_list$micro$SC$yearmo,
    df = ts_list$micro$SC,
    color = size_cols['micro']
  )+
  geom_point(
    data = post_chl$data$all_data |> 
      filter(frac == 'micro', site == 'SC'),
    aes(
      x = yearmo,
      y = ugL_chl
    )
  )+
  theme_pubclean()


# endregion -----------------------

# region Marginal Effect Plots -----------------------

marg_list <- list()


# region \- Add seasonal -----------------------

marg_list[['month']] <- list()
for(frac in size_factors) {
  marg_list[['month']][[frac]] <- list()
  frac_idx <- as.numeric(size_factors[which(size_factors == frac)])
  for(site in site_factors) {
    site_idx <- as.numeric(site_factors[which(site_factors == site)])
    sin_beta <- post_chl$beta[,frac_idx, site_idx,2]
    cos_beta <- post_chl$beta[,frac_idx, site_idx,3]
    amp <- sqrt(sin_beta^2 + cos_beta^2)
    phase <- atan2(cos_beta, sin_beta)

    chl_pred <- sapply(1:12, function(t) {
      post_chl$beta[,frac_idx, site_idx,1] +
        amp * sin(2*pi*t/12 + phase)
    }) |> 
      summarize_pred()

    chl_pred$month <- 1:12
    marg_list[['month']][[frac]][[site]] <- chl_pred
  }
  marg_list[['month']][[frac]] <- marg_list[['month']][[frac]] |> 
    EcotaxaTools::list_to_tib('site')
}

marg_list[['month']] <- marg_list[['month']] |> 
  EcotaxaTools::list_to_tib('frac')

# endregion -----------------------

# region Core Loop -----------------------

# loop for most variables
for(var in names(post_chl$data$X_scaled)[-c(1:3)]) {
  marg_list[[var]] <- list()
  var_idx <- which(names(post_chl$data$X_scaled) == var)
  for(frac in size_factors) {
    marg_list[[var]][[frac]] <- list()
    frac_idx <- as.numeric(size_factors[which(size_factors == frac)])
    for(site in site_factors) {
      site_idx <- as.numeric(site_factors[which(site_factors == site)])
      sub_data <- post_chl$data$all_data[which(
        post_chl$data$all_data$site == site & post_chl$data$all_data$frac == frac
      ),]
      sim_range <- seq(
        min(sub_data[[var]]), 
        max(sub_data[[var]]),
        length.out = 100
      )
      if(var %in% c('N',"NH4","P")) {
        sim_scale <- scale(
          log(sim_range +1e-5), 
          center = mean(log(post_chl$data$all_data[[var]]+1e-5)), 
          scale = sd(log(post_chl$data$all_data[[var]]+1e-5)))
      } else {
        sim_scale <- scale(
          sim_range,
          center = mean(post_chl$data$all_data[[var]]),
          scale = sd(post_chl$data$all_data[[var]])
        )
      }

      fit_chl <- post_chl$beta[,frac_idx, site_idx,1] +
        as.matrix(sim_scale) %*% post_chl$beta[,frac_idx, site_idx,var_idx]
      fit_chl <- summarize_pred(t(fit_chl))
      fit_chl[[var]] <- sim_range

      marg_list[[var]][[frac]][[site]] <- fit_chl
    }
    marg_list[[var]][[frac]] <- marg_list[[var]][[frac]] |> 
      EcotaxaTools::list_to_tib('site')
  }
  marg_list[[var]] <- marg_list[[var]] |> 
    EcotaxaTools::list_to_tib('frac')
}

# endregion -----------------------

ggplot(
  marg_list[[var]] |> 
    mutate(
      site = factor(site, levels = site_factors)
    )
) +
  geom_ribbon(
    aes(
      x =  marg_list[[var]][[var]],
      ymin = low.95,
      ymax = high.95,
      fill = frac
    ),
    alpha = 0.5
  ) +
  geom_ribbon(
    aes(
      x =  marg_list[[var]][[var]],
      ymin = low.75,
      ymax = high.75,
      fill = frac
    ),
    alpha = 0.5
  ) +
  geom_ribbon(
    aes(
      x =  marg_list[[var]][[var]],
      ymin = low.50,
      ymax = high.50,
      fill = frac
    ),
    alpha = 0.5
  ) +
 
  facet_grid(
    site ~ frac
  ) +
  scale_fill_manual(
    values = size_cols
  )+
  guides(fill = 'none')+
  labs(x = var, y = expression(paste("Chl-a ", "[",mu * g~L^{-1},"]")))+
  theme_pubclean() +
  theme(strip.background = element_rect(fill = 'transparent'))

# endregion -----------------------
