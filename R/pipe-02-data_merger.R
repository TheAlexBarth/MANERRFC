#########################
# Data merger -----------
#########################
rm(list = ls())
library(ggplot2)
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
    across(c(t, t_max, t_min, s_max, sal, DO, DO_min, Chl, Chl_max, Chl_min, Turb), \(x) mean(x, na.rm = T))
  )

# average to the yearmo for nuts
nut_yearmo <- environ$nut_avg |> 
  mutate(
    yearmo = paste(year(.data$date),month(.data$date),'01',sep = '-') |> 
      as.Date()
  ) |> 
  group_by(yearmo, sample_site) |> 
  summarize(
    across(c(P,NH4,N, CHLA_N), \(x) mean(x, na.rm = T))
  )

# average to the yearmo for winds
wind_yearmo <- environ$wind_avg |>
  mutate(
    yearmo = paste(year(.data$date),month(.data$date),'01',sep = '-') |>
      as.Date()
  ) |>
  group_by(yearmo) |>
  summarize(
    across(c(windspeed, TotalPAR), mean, na.rm = T)
  )

# region \- merge conc ----------


#attach variates to each concentration. Start with All

conc_mg <- etx$conc |> 
  left_join(
    wq_yearmo,
    by = c('sample_site', 'yearmo')
  ) |> 
  left_join(
    nut_yearmo,
    by = c('sample_site', 'yearmo')
  ) |>
  left_join(
    wind_yearmo,
    by = c('yearmo')
  )


# region \- merge indv. ----------------
indv_merged <- etx$indv |>
  left_join(
    wq_yearmo,
    by = c('sample_site', 'yearmo')
  ) |> 
  left_join(
    nut_yearmo,
      by = c('sample_site', 'yearmo')
  ) |> 
  left_join(
    wind_yearmo,
    by = c('yearmo')
  )

#  curiousity plot
#  ggplot(
#   conc_mg |> 
#     filter(taxa %in% etx$name$diatom)
#   ) +
#    geom_point(
#      aes(
#        x = windspeed,
#        y = log(pgC_L+1)
#      )
#    ) +
#    geom_abline(slope = 1, intercept = 0)


# # region \- merge living --------------
# living <- etx$indv |> 
#   left_join(
#     wq_yearmo,
#     by = c('sample_site', 'yearmo')
#   ) |> 
#   left_join(
#     nut_yearmo,
#     by = c('sample_site', 'yearmo')
#   ) |>
#    left_join(
#      wind_yearmo,
#      by = c('yearmo')
#    )


#  tot_conc <- conc_mg |> 
#    group_by(id,sample_site) |> 
#    summarize(
#      pgC_L = sum(pgC_L, na.rm = T)
#    ) |> 
#    left_join(
#      conc_mg |> 
#        select(id, t, Chl, sal),
#      by = 'id'
#    ) |> 
#    unique()


 
saveRDS(
  list(
    names = etx$names,
    conc = conc_mg,
    indv = indv_merged
  ),
  './data/02-full_merged.rds'
)
