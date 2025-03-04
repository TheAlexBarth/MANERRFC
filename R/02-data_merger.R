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

# average to the yearmo for nuts
nut_yearmo <- environ$nut_avg |> 
  mutate(
    yearmo = paste(year(.data$date),month(.data$date),'01',sep = '-') |> 
      as.Date()
  ) |> 
  group_by(yearmo, sample_site) |> 
  summarize(
    across(c(P,NH4,N, CHLA_N), mean, na.rm = T)
  )

# region \- merge conc ----------

conc_mg <- etx$conc |> 
  left_join(
    wq_yearmo,
    by = c('sample_site', 'yearmo')
  ) |> 
  left_join(
    nut_yearmo,
    by = c('sample_site', 'yearmo')
  )

# # curiousity plot
# ggplot(conc_mg) +
#   geom_point(
#     aes(
#       x = Chl,
#       y = CHLA_N,
#       color = sample_site
#     )
#   ) +
#   geom_abline(slope = 1, intercept = 0)


# region \- merge living --------------
living <- etx$indv |> 
  left_join(
    wq_yearmo,
    by = c('sample_site', 'yearmo')
  ) |> 
  left_join(
    nut_yearmo,
    by = c('sample_site', 'yearmo')
  )


# tot_conc <- conc_mg |> 
#   group_by(id,sample_site) |> 
#   summarize(
#     pgC_L = sum(pgC_L, na.rm = T)
#   ) |> 
#   left_join(
#     conc_mg |> 
#       select(id, t, Chl, sal),
#     by = 'id'
#   ) |> 
#   unique()

# ggplot(conc_mg) +
#   geom_point(
#     aes(
#       x = sal,
#       y = Chl,
#       color = sample_site
#     )
#   ) +
#   geom_smooth(
#     aes(
#       x = sal,
#       y = Chl,
#       color = sample_site
#     ),
#     method = 'lm'
#   )+
#   theme_minimal()

saveRDS(
  list(
    names = etx$names,
    conc = conc_mg,
    indv = living
  ),
  './data/02-full_merged.rds'
)