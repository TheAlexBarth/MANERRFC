rm(list = ls())
library(EcotaxaTools) # this is my package
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(rstan)
library(lubridate)
set.seed(0121)
source('./R/utils.R')


# read conc data
conc <- readRDS('./data/t01-copano_concentrations.rds')
name_list = readRDS('./data/t01-taxa_names.rds')

# read individual data
indv <- readRDS('./data/t01-copano_individals.rds')


# read wq
wq = readRDS('./output/s02-temp-wq.rds')


# region \- water quality prep for ts ------------------

wq = wq |> 
    filter(sample_site %in% c('CE', 'CW')) |> 
    group_by(Date) |> 
    summarise(sal = mean(Sal), chl = mean(ChlFluor, na.rm = T))

wq = wq[!is.na(wq$sal),]


########################
# MARK: Diatoms --------
########################


# trim to diatoms
diat_conc <- conc |> 
  filter(taxa %in% name_list$diatom)

# flip taxa names 
diat_conc$group <- diat_conc$taxa |> 
  sapply(
    function(x) names(name_list$diatom)[which(name_list$diatom == x)]
  ) |> 
  as.factor()


all_diat <- diat_conc |>
  group_by(sample_id, group) |> 
  summarize(count = sum(count), num_L = sum(num_L), um3_L = sum(um3_L)) |> 
  left_join(
    diat_conc |> 
    select(-c(count, num_L, um3,um3_L, taxa, group)) |> 
      unique(),
    by = "sample_id"
  )

# impute to mean
all_diat$sal[which(is.na(all_diat$sal))] <- mean(all_diat$sal, na.rm = T)

diat_indv <- indv |> 
  filter(taxo_name %in% name_list$diatom)

diat_indv$group <- diat_indv$taxo_name |> 
  sapply(
    function(x) names(name_list$diatom)[which(name_list$diatom == x)]
  ) |> 
  as.factor()


# region \- stan fit ----------
diat_indv$sal[is.na(diat_indv$sal)] <- mean(diat_indv$sal, na.rm = T)

# don't use mix-multiple since there is no effect of temperature
diat_pois <- list(
  N_obs = nrow(all_diat),
  N_mes = nrow(diat_indv),
  N_groups = length(unique(all_diat$group)),
  n = all_diat$count,
  group = as.numeric(all_diat$group),
  img_vol = all_diat$img_vol/1000, # convert to L,
  sal = all_diat$sal,
  group_counts = as.numeric(all_diat$group),
  group_bio = as.numeric(diat_indv$group),
  log_b = log(diat_indv$um3),
  sal_b = diat_indv$sal,
  temp_b = diat_indv$t,
  obs_biomass_conc = all_diat$um3_L
)

diat_fit <- stan(
  file = './stan/03b-mix-sal_only-hierachical.stan',
  data = diat_pois,
  chains = 4, iter = 3000, warmup = 500, cores = 14
)

# region \- prediction --------------------------

# predict lambda
g_ts_list <- list()
lamba_ts <- list()
n_ts <- list()
for(l in 1:length(levels(all_diat$group))) {
  lamba_ts[[l]] <- post_predict(
    extract(diat_fit),
    list(
      sal = wq$sal
    ),
    paste0(
      'exp(beta0[,',l,'] + beta1[,',l,'] * sal + log(0.1129333/1000))'
    ),
    n_draws = 2000
  )

  n_ts[[l]] <- matrix(
    rpois(length(lamba_ts[[l]]), lambda = lamba_ts[[l]]),
    nrow = nrow(lamba_ts[[l]]),
    ncol = ncol(lamba_ts[[l]])
  )
}

# Predict Biomass
err_ts <- extract(diat_fit, pars = 'v')$v
eta_ts <- list()
d_ts <- list()
for(l in 1:length(levels(all_diat$group))) {
  eta_ts[[l]] <- post_predict(
    extract(diat_fit),
    list(
      sal = wq$sal
    ),
    paste0(
      'a0[,',l,'] + a1[,',l,'] * sal'
    ),
    n_draws = 2000
  )
 
  d_ts[[l]] <- apply(
    eta_ts[[l]], 2, 
    function(x) rnorm(nrow(eta_ts[[l]]), mean = x, sd = err_ts)
  ) |> 
    exp()

  g_ts_list[[l]] = d_ts[[l]] * (n_ts[[l]]/mean(all_diat$img_vol/1000))
}

# region \- total fits ---------------


total_diat <- all_diat |> 
  group_by(sample_id, date) |> 
  summarize(
    biomass = sum(um3_L, na.rm = T),
    sal = mean(sal, na.rm = T)
  )


total_ts <- g_ts_list |> 
  summarize_list(sp_range = wq$sal)

total_plot_ts <- wq |> 
  left_join(
    total_ts,
    by = 'sal'
  )


ggplot() +
  geom_point(
    data = total_diat,
    aes(
      x = date,
      y = biomass
    ),
    size = 3, color = '#f2dbf6'
  )+
  geom_error_range(total_plot_ts$Date, total_plot_ts, '#cc91ff') +
  scale_x_date(
    limits = c(min(total_plot_ts$Date), max(total_plot_ts$Date))
  )+
  labs(x = "", y = '')+
  scale_y_continuous(limits = c(0,6.5e10))+
  theme_minimal()+
  theme(
    axis.text = element_text(size = 28, color = 'white'),
    panel.grid = element_line(color = 'grey50'),
    plot.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.background = element_rect(fill = 'transparent')
  )

  ggsave('./output/msac/diat-ts.pdf',
  height = 6, width = 22, units = 'in')
        