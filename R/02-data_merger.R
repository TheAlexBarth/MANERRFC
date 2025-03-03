#########################
# Data merger -----------
#########################
rm(list = ls())
library(dplyr)
library(lubridate)

etx = readRDS('./data/00-ecotaxa_full.rds')
environ = readRDS('./data/01-environ_clean.rds')


# water quality should be summed in month
wq_yearmo <- environ$wq_sum |> 
  mutate(
    yearmo = paste(year(.data$date),month(.data$date),'01',sep = '-') |> 
      as.Date()
  ) |> 
  group_by(yearmo, sample_site) |> 
  summarize(
    across(c(t, t_max, t_min, s_max, sal, DO, DO_min, Chl, Chl_max, Chl_min), mean, na.rm = T)
  )


# region \- for conc ----------

