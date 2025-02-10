rm(list = ls())
library(EcotaxaTools) # this is my package
library(dplyr)
library(ggplot2)
library(tidyr)
library(rstan)
library(lubridate)
library(ggpubr)
source('./R/utils.R')

## if time need to run as non-hierarchical
conc <- readRDS('./data/t01-copano_concentrations.rds')
name_list = readRDS('./data/t01-taxa_names.rds')
indv <- readRDS('./data/t01-copano_individals.rds')

# trim to dinos
conc <- conc |> 
  filter(taxa %in% name_list$dino)

# flip taxa names 
conc$group <- conc$taxa |> 
  sapply(
    function(x) names(name_list$dino)[which(name_list$dino == x)]
  ) |> 
  as.factor()



#####
# MARK: data prep 
######

all_dino <- conc |>
  group_by(sample_id, group) |> 
  summarize(count = sum(count), num_L = sum(num_L), um3_L = sum(um3_L)) |> 
  left_join(
    conc |> 
    select(-c(count, num_L, um3,um3_L, taxa, group)) |> 
      unique(),
    by = "sample_id"
  )

# impute to mean
all_dino$sal[which(is.na(all_dino$sal))] <- mean(all_dino$sal, na.rm = T)
all_dino$t[which(is.na(all_dino$t))] <- mean(all_dino$t, na.rm = T)
# individual metrics
indv_metrics <- indv |> 
  filter(taxo_name %in% name_list$dino)

indv_metrics$group <- indv_metrics$taxo_name |> 
  sapply(
    function(x) names(name_list$dino)[which(name_list$dino == x)]
  ) |> 
  as.factor()
indv_metrics$sal[is.na(indv_metrics$sal)] <- mean(indv_metrics$sal, na.rm = T)
indv_metrics$t[is.na(indv_metrics$t)] <- mean(indv_metrics$t, na.rm = T)




# region \- quick plots

ggplot() +
  geom_point(
    data = all_dino,
    aes(
      x = sal,
      y = um3_L,
      color = group
    )
  )


####
# MARK: MODEL 
####


# region \- stan fit ----------

# don't use mix-multiple since there is no effect of temperature
pois_data <- list(
  N_obs = nrow(all_dino),
  N_mes = nrow(indv_metrics),
  N_groups = length(unique(all_dino$group)),
  n = all_dino$count,
  group = as.numeric(all_dino$group),
  img_vol = all_dino$img_vol/1000, # convert to L,
  sal = all_dino$sal,
  temp = all_dino$t,
  group_counts = as.numeric(all_dino$group),
  group_bio = as.numeric(indv_metrics$group),
  log_b = log(indv_metrics$um3),
  sal_b = indv_metrics$sal,
  temp_b = indv_metrics$t,
  obs_biomass_conc = all_dino$um3_L
)

pois_fit <- stan(
  file = './stan/03c-mix-multiple.stan',
  data = pois_data,
  chains = 4, iter = 3000, warmup = 500, cores = 14
)

#
####################
# MARK: Single post evals
####################
sal_pred = seq(min(all_dino$sal, na.rm = T), max(all_dino$sal, na.rm = T), by = 0.1)
temp_pred = seq(min(all_dino$t, na.rm = T), max(all_dino$t, na.rm = T), by = 0.1)

# region: \- sal
sal_pred_list <- list()
for(l in 1:length(levels(all_dino$group))) {
  sal_pred_list[[l]] <- post_predict(
    extract(pois_fit),
    list(
      sal = sal_pred,
      temp = rep(mean(temp_pred), length(sal_pred))
    ),
    paste0(
      'exp(beta0[,',l,'] + beta1[,',l,'] * sal + beta2[,',l,'] * temp) * exp(a0[,',l,'])'
    ),
    n_draws = 2000
  )
}

marg <- sal_pred_list[[6]] |> 
  summarize_pred()

ggplot() +
  # geom_point(
  #   data = all_dino |> 
  #     filter(group == 'Margalef'),
  #   aes(
  #     x = sal,
  #     y = um3_L
  #   )
  # ) +
  geom_line(
    aes(
      x = sal_pred,
      y = marg$mean
    )
  ) +
  geom_ribbon(
    aes(
      x = sal_pred,
      ymin = marg$low,
      ymax = marg$high
    ),
    alpha = 0.25
  )

# region: \- temp
t_pred_list <- list()
for(l in 1:length(levels(all_dino$group))) {
  t_pred_list[[l]] <- post_predict(
    extract(pois_fit),
    list(
      sal = rep(mean(all_dino$sal), length(temp_pred)),
      temp = temp_pred
    ),
    paste0(
      'exp(beta0[,',l,'] + beta1[,',l,'] * sal + beta2[,',l,'] * temp) * exp(a0[,',l,'])'
    ),
    n_draws = 2000
  )
}



# region \- marg
marg <- t_pred_list[[1]] |> 
  summarize_pred()

ggplot() +
  geom_point(
    data = all_dino |> 
      filter(group == 'Akashiwo'),
    aes(
      x = t,
      y = um3_L
    )
  ) +
  geom_line(
    aes(
      x = temp_pred,
      y = marg$mean
    )
  ) +
  geom_ribbon(
    aes(
      x = temp_pred,
      ymin = marg$low,
      ymax = marg$high
    ),
    alpha = 0.25
  )


############
# MARK: MEAN TREND 
##############


# region: \- sal
mean_trend <- post_predict(
  extract(pois_fit),
  list(
    sal = sal_pred,
    temp = rep(mean(temp_pred), length(sal_pred))
  ),
  paste0(
    'exp(mu0+ mu1 * sal + mu2 * temp) * exp(theta0)'
  ),
  n_draws = 3000
) |> 
  summarize_pred()

ggplot() +
  geom_point(
    data = all_dino,
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
  scale_color_manual(values = gg_cbb_col(length(levels(all_dino$group))))+
  theme_pubclean() +
  theme(legend.position = 'none')

# region: \- temp

mean_trend <- post_predict(
  extract(pois_fit),
  list(
    sal = rep(mean(all_dino$sal), length(temp_pred)),
    temp = temp_pred
  ),
  paste0(
    'exp(mu0 + mu1 * sal + mu2 * temp) * exp(theta0)'
  ),
  n_draws = 300
) |> 
  summarize_pred()

ggplot() +
  geom_point(
    data = all_dino,
    aes(
      x = t,
      y = um3_L,
      color = group
    )
  ) +
  geom_line(
    aes(
      x = temp_pred,
      y = mean_trend$mean
    )
  ) +
  geom_ribbon(
    aes(
      x = temp_pred,
      ymin = mean_trend$low,
      ymax = mean_trend$high
    ),
    alpha = 0.25
  ) +
  scale_color_manual(values = gg_cbb_col(length(levels(all_dino$group))))+
  theme_pubclean() +
  theme(legend.position = 'none')


# region \- summarize group
sum_pred <- t_pred_list |> 
  lapply(summarize_pred) |> 
  lapply(
    function(x) mutate(x, t = temp_pred)
  )
names(sum_pred) <- levels(all_dino$group)
total_pred <- sum_pred |> 
  list_to_tib() |> 
  group_by(t) |> 
  summarize(
    mean = sum(mean),
    low = sum(low),
    high = sum(high)
  )

total_dino <- all_dino |> 
  group_by(sample_id) |> 
  summarize(
    biomass = sum(um3_L, na.rm = T),
    temp = mean(t, na.rm = T),
    sal = mean(sal, na.rm = T)
  )

ggplot() +
  geom_point(
    aes(
      x = total_dino$sal,
      y = total_dino$biomass
    )
  )
  geom_line(
    data = total_pred,
    aes(
      x = t,
      y = mean
    )
  ) +
  geom_ribbon(
    data = total_pred,
    aes(
      x = t,
      ymin = low,
      ymax = high
    ),
    alpha = 0.25
  ) +
  scale_y_log10()
