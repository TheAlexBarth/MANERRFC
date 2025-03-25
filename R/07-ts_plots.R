rm(list = ls())
library(ggplot2)
library(tidyr)
library(dplyr)

# you need to run sequentially
# set the taxa in first chunk
# interactively make plots at bottom

source('./R/utils.R')
source('./R/utils-posterior_specific.R')
source('./R/utils-mod_making.R')

###########
# MARK: CHOOSE TAXA 
#############

taxa = 'diatom'


# THIS ONLY WORKS BECAUSE NONE ARE INFLUENCED BY NUTRIENTS

##############################
# MARK: LOAD DATA
##############################


etx <- readRDS('./data/02-full_merged.rds')
group_data <- prep_component_data(taxa)
mod = readRDS(paste0('./data/04-mod_output-',taxa,'.RDS'))
wq_full = readRDS('./data/01-environ_clean.rds')

wq <- wq_full$wq_sum |> 
  select(date, sample_site, t, sal, do = DO, turb = Turb) |> 
  left_join(
    wq_full$wind_avg |> 
      select(date, windspeed),
    by = 'date'
  )



########################
# MARK: Posterior Format
########################



post_coef = rstan::extract(mod$mod, pars = c('beta_bio','beta_count', 'v'))

########################
# MARK: prepare data
########################

# need to do for each site

ce_pred <- make_ts_data('CE', wq = wq, data = group_data, mod = mod, post_coef = post_coef)


#################
# MARK: MAKE PLOT
#################


ce_ts_total <- make_ts_plot(
  taxa = 'Thalassionema',
  site = 'CE',
  sim_data = ce_pred,
  fill_col = '#4500f4',
  size = 2 #affect the points
)

ce_ts_total +
  labs(x = "", y = "pgC / L") +
  scale_x_date()+
  theme_minimal()
