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

###########################
# MARK: data prep --------
###########################

# read conc data
conc <- readRDS('./data/t01-copano_concentrations.rds')
name_list = readRDS('./data/t01-taxa_names.rds')

# read individual data
indv <- readRDS('./data/t01-copano_individals.rds')


# read wq
wq = readRDS('./output/s02-temp-wq.rds')
nut <- readRDS('./data/01_swmp_nut.rds')

# region \- water quality prep for ts ------------------

wq = wq |> 
    filter(sample_site %in% c('CE', 'CW')) |> 
    group_by(Date) |> 
    summarise(
      sal = mean(Sal), chl = mean(ChlFluor, na.rm = T),
      temp = mean(temp)
    )

wq = wq[!is.na(wq$sal),]



########################
# MARK: Chlorophyll ---
########################

# region \-data prep --------------
nut_copano <- nut |> 
  filter(substr(nut$StationCode, 4,5) %in% c('ce','cw'))

nut_copano$Date <- as.Date(nut_copano$DateTimeStamp)

nut_full <- nut_copano |> 
  select(Date, chl = CHLA_N) |> 
  group_by(Date) |> 
  summarize(chl = mean(chl)) |> 
  left_join(wq |> select(Date, sal), by = 'Date') |> 
  filter(!is.na(chl) & !is.na(sal))

# region \-stan ----


chl_stan <- list(
  N = nrow(nut_full),
  K = 2,
  x_covars = cbind(1, nut_full$sal),
  y = nut_full$chl
)

chl_mod <- stan(
  './stan/t0121-simple_reg.stan',
  data = chl_stan
)

# region \- post pred----------------------------
chl_beta <- extract(chl_mod, pars = 'beta')$beta |> 
  t()

sal_range <- seq(0,42, 0.1)

pred_chl_mean <- cbind(1, sal_range) %*% chl_beta
pred_chl_sigma <- extract(chl_mod, pars = 'sigma')$sigma |> mean()


pred_chl <- rnorm(
  length(pred_chl_mean), mean = pred_chl_mean, sd = pred_chl_sigma
) |> 
  matrix(nrow = nrow(pred_chl_mean), ncol = ncol(pred_chl_mean))

pred_chl <- pred_chl |> t() |> summarize_pred()

plot_mean <- pred_chl_mean |> t() |> summarize_pred()

 # region \- plot ---------------------
ggplot() +
  geom_point(
    data = nut_full,
    aes(
      x = sal,
      y = chl
    ),
    color = '#49be31'
  ) +
  geom_error_range(
    x = sal_range,
    plot_mean,
    color = '#90f87b'
  )+
  scale_y_continuous(limits = c(0,25))+
  scale_x_continuous(limits = c(0,45))+
  labs(x = "", y = '')+
  theme_minimal() +
  theme(
    axis.text = element_text(size = 28, color = 'white'),
    panel.grid = element_line(color = 'grey50'),
    plot.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.background = element_rect(fill = 'transparent')
  )

ggsave('./output/msac/01-chl.pdf',
height = 8, width = 8, units = 'in')







########################
# MARK: Time-series ---
########################


# region \- pred chl -------------------
ts_chl_mean <- cbind(1, wq$sal) %*% chl_beta

pred_chl <- rnorm(
  length(ts_chl_mean), mean = ts_chl_mean, sd = pred_chl_sigma
) |> 
  matrix(nrow = nrow(ts_chl_mean), ncol = ncol(ts_chl_mean))

pred_chl <- pred_chl |> t() |> summarize_pred()


# region \- plot ------------
ggplot() +
  geom_point(
    aes(
      x = nut_full$Date,
      y = nut_full$chl
    ),
    color = '#49be31'
  )+
  geom_error_range(
    x = wq$Date,
    pred_chl,
    '#90f87b'
  )+
  scale_x_date(
    limits = c(min(wq$Date), max(wq$Date))
  )+
  labs(x = "", y = '')+
  theme_minimal()+
  theme(
    axis.text = element_text(size = 28, color = 'white'),
    panel.grid = element_line(color = 'grey50'),
    plot.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.background = element_rect(fill = 'transparent')
  )


ggsave('./output/msac/chl-ts.pdf',
height = 6, width = 22, units = 'in')
      
