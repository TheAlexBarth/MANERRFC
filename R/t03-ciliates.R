rm(list = ls())
library(EcotaxaTools) # this is my package
library(dplyr)
library(ggplot2)
library(tidyr)
library(rstan)
library(brms)
library(tidybayes)
library(emmeans)
library(ggpubr)
library(lubridate)
source('./R/utils.R')


conc <- readRDS('./data/t01-copano_concentrations.rds')
name_list = readRDS('./data/t01-taxa_names.rds')
indv <- readRDS('./data/t01-copano_individals.rds')



##################
# MARK: mzom Total
##################

## \- mzoms --------

all_mz <- conc |>
  filter(taxa %in% name_list$mz) |> 
  group_by(sample_id) |> 
  summarize(count = sum(count), num_L = sum(num_L), um3_L = sum(um3_L)) |> 
  left_join(
    conc |> 
    select(-c(count, num_L, um3,um3_L, taxa)) |> 
      unique(),
    by = "sample_id"
  )

# impute to mean
all_mz$sal[which(is.na(all_mz$sal))] <- mean(all_mz$sal, na.rm = T)
# for later plots
sal_pred = seq(min(all_mz$sal, na.rm = T), max(all_mz$sal, na.rm = T), by = 0.1)


# region \- simple plots

plot(all_mz$um3_L ~ all_mz$sal)
abline(lm(all_mz$um3_L ~ all_mz$sal))

plot(log10(all_mz$um3_L+1) ~ all_mz$sal)
abline(lm(log10(all_mz$um3_L+1) ~ all_mz$sal))


indv_metrics <- indv |> 
  filter(taxo_name %in% name_list$diatom)

plot(log(all_mz$count * mean(indv_metrics$um3)) ~ all_mz$sal)

ggplot() +
  geom_point(
    aes(
      x = all_mz$sal,
      y = (all_mz$count * mean(indv_metrics$um3))
    )
  )

ggplot() +
  geom_point(
    aes(
      x = all_mz$sal,
      y = (all_mz$count * mean(indv_metrics$um3))
    )
  )


# region \- stan fit ----------
indv_metrics$sal[is.na(indv_metrics$sal)] <- mean(indv_metrics$sal, na.rm = T)
pois_data <- list(
  N_obs = nrow(all_mz),
  N_mes = nrow(indv_metrics),
  n = all_mz$count,
  img_vol = all_mz$img_vol/1000, # convert to L,
  sal = all_mz$sal,
  log_b = log(indv_metrics$um3),
  sal_b = indv_metrics$sal,
  temp_b = indv_metrics$t,
  obs_biomass_conc = all_mz$um3_L
)

mult_fit <- stan(
  file = './stan/03-mix-sal_only.stan',
  data = pois_data,
  chains = 4, iter = 4000, warmup = 500, cores = 14
)

pois_draws <- extract(mult_fit, pars = c('beta0', 'beta1','a0'))
pred_counts <- post_predict(
  pois_draws, 
  list(sal = sal_pred), 
  "exp(beta0 + beta1 * sal)"
) |> 
  summarize_pred()


# region \- Plotting


ggplot() +
  geom_point(
    aes(
      x = all_mz$sal,
      y = all_mz$num_L
    ),
    size = 2.3,
    alpha = 0.5
  ) +
  geom_line(
    aes(
      x = sal_pred,
      y = pred_counts$mean
    ),
    color = 'red',
    linewidth = 2
  ) +
  geom_ribbon(
    aes(
      x = sal_pred,
      ymin = pred_counts$low,
      ymax = pred_counts$high
    ),
    fill = 'red',
    alpha = 0.25
  ) +
  labs(x = "Salinity", y = "Num / L")+
  theme_bw()

ggsave('./output/t03-ciliate-count.pdf', width = 3, height = 3, units = 'in', dpi = 600)



pred_biov <-  post_predict(
  pois_draws, 
  list(sal = sal_pred), 
  "exp(beta0 + beta1 * sal) * exp(a0)"
) |> 
  summarize_pred()



ggplot() +
  geom_point(
    aes(
      x = all_mz$sal,
      y = all_mz$um3_L +1
    )
  )   +
  geom_line(
    aes(
      x = sal_pred,
      y = pred_biov$mean
    ),
    color = 'red'
  ) +
  geom_ribbon(
    aes(
      x = sal_pred,
      ymin = pred_biov$low,
      ymax = pred_biov$high
    ),
    fill = 'red',
    alpha = 0.25
  ) +
  labs(x = "Salinity", y = "um3 / L")+
  theme_bw()

##############################
# MARK: Multiple Groups
##############################
names(name_list$mz)[which(names(name_list$mz) == 'dead_choreo')] = 'Tintinnina'


g_conc <- conc |> 
  filter(taxa %in% name_list$mz)

g_conc$group <- g_conc$taxa |> 
  sapply(
    function(x) names(name_list$mz)[which(name_list$mz == x)]
  ) |> 
  as.factor()

g_mz <- g_conc |>
  group_by(sample_id, group) |> 
  summarize(count = sum(count), num_L = sum(num_L), um3_L = sum(um3_L)) |> 
  left_join(
    g_conc |> 
    select(-c(count, num_L, um3,um3_L, taxa, group)) |> 
      unique(),
    by = "sample_id"
  )

# impute
g_mz$sal[which(is.na(g_mz$sal))] <- mean(g_mz$sal, na.rm = T)
g_mz$t[which(is.na(g_mz$t))] <- mean(g_mz$t, na.rm = T)


indv_metrics <- indv |> 
  filter(taxo_name %in% name_list$mz)

indv_metrics$group <- indv_metrics$taxo_name |> 
  sapply(
    function(x) names(name_list$mz)[which(name_list$mz == x)]
  ) |> 
  as.factor()
indv_metrics$sal[is.na(indv_metrics$sal)] <- mean(indv_metrics$sal, na.rm = T)
indv_metrics$t[is.na(indv_metrics$t)] <- mean(indv_metrics$t, na.rm = T)

# region \- quick plots

ggplot() +
  geom_point(
    data = g_mz,
    aes(
      x = sal,
      y = um3_L,
      color = group
    )
  )

ggplot(    data = indv_metrics,
  aes(
    x = t,
    y = log(um3)
  )) +
  geom_point(

  ) +
  geom_smooth(method = 'lm')

ggplot(    data = g_mz,
  aes(
    x = sal,
    y = t
  )) +
  geom_point(

  ) +
  geom_smooth(method = 'lm')

# region \- stan fit
# don't use mix-multiple since there is no effect of temperature
mult_data <- list(
  N_obs = nrow(g_mz),
  N_mes = nrow(indv_metrics),
  N_groups = length(unique(g_mz$group)),
  n = g_mz$count,
  group = as.numeric(g_mz$group),
  img_vol = g_mz$img_vol/1000, # convert to L,
  sal = g_mz$sal,
  group_counts = as.numeric(g_mz$group),
  group_bio = as.numeric(indv_metrics$group),
  log_b = log(indv_metrics$um3),
  sal_b = indv_metrics$sal,
  temp_b = indv_metrics$t,
  obs_biomass_conc = g_mz$um3_L
)

mult_fit <- stan(
  file = './stan/03b-mix-sal_only-hierachical.stan',
  data = mult_data,
  chains = 4, iter = 3000, warmup = 500, cores = 14
)


########################
# MARK: POSTERIOR EVALS
########################
sal_pred = seq(min(g_mz$sal, na.rm = T), max(g_mz$sal, na.rm = T), by = 0.1)


# predict lambda
g_pred_list <- list()
lambda_fit <- list()
n_pred <- list()
for(l in 1:length(levels(g_mz$group))) {
  lambda_fit[[l]] <- post_predict(
    extract(mult_fit),
    list(
      sal = sal_pred
    ),
    paste0(
      'exp(beta0[,',l,'] + beta1[,',l,'] * sal + log(0.1129333/1000))'
    ),
    n_draws = 2000
  )

  n_pred[[l]] <- matrix(
    rpois(length(lambda_fit[[l]]), lambda = lambda_fit[[l]]),
    nrow = nrow(lambda_fit[[l]]),
    ncol = ncol(lambda_fit[[l]])
  )
}

# Predict Biomass
err_var <- extract(mult_fit, pars = 'v')$v
eta_pred <- list()
d_pred <- list()
for(l in 1:length(levels(g_mz$group))) {
  eta_pred[[l]] <- post_predict(
    extract(mult_fit),
    list(
      sal = sal_pred
    ),
    paste0(
      'a0[,',l,'] + a1[,',l,'] * sal'
    ),
    n_draws = 2000
  )
 
  d_pred[[l]] <- apply(
    eta_pred[[l]], 2, 
    function(x) rnorm(nrow(eta_pred[[l]]), mean = x, sd = err_var)
  ) |> 
    exp()

  g_pred_list[[l]] = d_pred[[l]] * (n_pred[[l]]/mean(g_mz$img_vol/1000))
}


# region \- plotting list specific
tin_plot <- g_pred_list[[7]] |> 
  summarize_pred()

ggplot() +
  geom_point(
    data = g_mz |> 
      filter(group == 'Tintinnina'),
    aes(
      x = sal,
      y = um3_L
    ),
    color = gg_cbb_col(2)[2],
    size = 3
  ) +
  geom_line(
    data = tin_plot,
    aes(
      x = sal_pred,
      y = mean
    ),
    color = gg_cbb_col(2)[2],
    linewidth = 1.5
  )+  
  geom_ribbon(
    data = tin_plot,
    aes(
      x = sal_pred,
      ymin = low,
      ymax = high
    ),
    alpha = 0.25,
    fill = gg_cbb_col(2)[2]
  )+
  labs(x = "Salinity", y = "Biovolume Concentration [um3/L]")+
  theme_pubclean()+
  theme(legend.position = 'none') +
  theme(
    axis.title = element_text(face = 'bold', size = 16),
    axis.text = element_text(face = 'bold', size = 12)
  )
      
    
    



# region \- mean trend
mean_trend <- post_predict(
  extract(mult_fit),
  list(
    sal = sal_pred
  ),
  paste0(
    'exp(mu0 + mu1 * sal) * exp(theta0 + theta1 * sal)'
  ),
  n_draws = 2000
) |> 
  summarize_pred()

ggplot() +
  geom_point(
    data = g_mz,
    aes(
      x = sal,
      y = um3_L,
      color = group
    )
  ) +
  geom_line(
    aes(
      x = sal_pred,
      y = mean_trend$mean
    )
  ) +
  geom_ribbon(
    aes(
      x = sal_pred,
      ymin = mean_trend$low,
      ymax = mean_trend$high
    ),
    alpha = 0.25
  ) +
  scale_color_manual(values = gg_cbb_col(length(levels(g_mz$group))))+
  theme_pubclean()+
  theme(legend.position = 'none')



# region \- summarize group
sum_pred <- g_pred_list |> 
  lapply(summarize_pred) |> 
  lapply(
    function(x) mutate(x, sal = sal_pred)
  )
names(sum_pred) <- levels(g_mz$group)
total_pred <- sum_pred |> 
  list_to_tib() |> 
  group_by(sal) |> 
  summarize(
    mean = sum(mean),
    low = sum(low),
    high = sum(high)
  )

total_mz <- g_mz |> 
  group_by(sample_id) |> 
  summarize(
    biomass = sum(um3_L, na.rm = T),
    sal = mean(sal, na.rm = T)
  )

ggplot() +
  geom_point(
    aes(
      x = total_mz$sal,
      y = total_mz$biomass
    ),
    color = "#b2b68c",
    size = 3
  )+
  geom_line(
    data = total_pred,
    aes(
      x = sal,
      y = mean
    ),
    color = "#b2b68c",
    linewidth = 1.5
  ) +
  geom_ribbon(
    data = total_pred,
    aes(
      x = sal,
      ymin = low,
      ymax = high
    ),
    fill = "#b2b68c",
    alpha = 0.25
  ) +
  labs(x = "Salinity", y = "Biovolume Concentration [um3/L]")+
  theme_pubclean()+
  theme(legend.position = 'none') +
  theme(
    axis.title = element_text(face = 'bold', size = 16),
    axis.text = element_text(face = 'bold', size = 12)
  )
  


# MARK: PLOT FOR SIZE #######

size_pred <- list()

for(l in 1:length(levels(g_mz$group))) {
  eta_pred[[l]] <- post_predict(
    extract(mult_fit),
    list(
      sal = 0
    ),
    paste0(
      'a0[,',l,']'
    ),
    n_draws = 2000
  )
 
  size_pred[[l]] <- apply(
    eta_pred[[l]], 2, 
    function(x) rnorm(nrow(eta_pred[[l]]), mean = x, sd = err_var)
  ) |> 
    exp() 

  size_pred[[l]] <- data.frame(
    biomass = size_pred[[l]][,1]
  )
}

names(size_pred) <- levels(g_mz$group)
size_pred <- size_pred |> 
  list_to_tib()

ggplot() +
  geom_density(
    data = size_pred,
    aes(
      x = log(biomass),
      fill = group
    ),
    alpha = 0.5
  ) +
  scale_fill_manual(values = rev(gg_cbb_col(8)[2:8])) +
  theme_pubclean() +
  labs(x = "Log(um3)", y = "", fill = "Group")
