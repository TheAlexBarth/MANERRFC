#########################
# Data merger -----------
#########################
rm(list = ls())
library(ggplot2)
library(dplyr)
library(lubridate)
source('./R/utils.R')

etx = readRDS('./data/00-ecotaxa_full.rds')

wq <- readRDS('./data/01a-swmp_wq_data.rds')
chl <- readRDS('./data/01b-size_frac_chl.RDS')
wind <- readRDS('./data/01c-wind_score.RDS')
si <- readRDS('./data/01e-silicate.rds')

# region \- merge conc ----------

# temp convert trophic role

#attach variates to each concentration. Start with All

troph_conc <- etx$conc |>
  group_by(functional_role, sample_site, yearmo, id) |>
  summarize(
    count = sum(count),
    img_vol = unique(img_vol)
  ) |>
  ungroup() |>
  left_join(
    wq |>
      select(P, NH4, N, CHLA_N, temp, sal, sampling_site, yearmo),
    by = c('sample_site' = 'sampling_site','yearmo')
  ) |>
  left_join(
    chl |>
      select(micro, nano, pico, site, yearmo),
     by = c('sample_site' = 'site','yearmo')
  ) |>
  left_join(
    wind |>
      select(wind_pca, yearmo),
     by = c('yearmo'),
     relationship = 'many-to-many'
  ) |>
  left_join(
    si |>
      select(SiOH, sampling_site, yearmo),
    by = c('sample_site' = 'sampling_site','yearmo')
  )

troph_conc$sample_site <- factor(troph_conc$sample_site, levels = levels(site_factors))
troph_conc$functional_role <- factor(troph_conc$functional_role, levels = levels(troph_factors))


# region \- merge indv. ----------------
indv_clean <- etx$indv |>
  select(
    yearmo, functional_role, um3, cmass, group, taxo_name, sample_site
  )

indv_clean$functional_role[which(grepl('grazer',indv_clean$functional_role))] <- 'grazer'

saveRDS(
  list(
    conc = troph_conc,
    indv = indv_clean
  ),
  './data/02-full_merged.rds'
)
