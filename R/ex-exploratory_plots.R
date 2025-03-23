########
# Exploratory Plots
#######

rm(list = ls())
library(ggplot2)
library(dplyr)
library(lubridate)
library(EcotaxaTools)


etx <- readRDS('./data/02-full_merged.rds')

###################################
# MARK: TIME SERIES STUFF -----
###################################

# region \- full time series plot ---------------------
# I make functions for plots that will be frequently reproduced
# but that's not required
full_ts_plotter <- function(data) {
  p = ggplot(data) +
    geom_point(
      aes(
        x = date,
        y = pgC_L,
        color = sample_site
      )
    ) +
    theme_minimal()
  return(p)
}

etx$conc |> 
  filter(taxa %in% etx$names$diatom) |> 
  group_by(date, sample_site) |> 
  summarize(pgC_L = sum(pgC_L)) |> 
  full_ts_plotter()

etx$conc |> 
  filter(taxa %in% etx$names$mz) |> 
  group_by(date, sample_site) |> 
  summarize(pgC_L = sum(pgC_L)) |> 
  full_ts_plotter()



etx$conc |> 
  filter(taxa == 'Spirotrichea') |> 
  group_by(date, sample_site) |> 
  summarize(pgC_L = sum(pgC_L)) |> 
  full_ts_plotter()


# region \- seasonality ------------
seasonality_plot <- function(groups) {
  data = etx$conc |> 
    filter(taxa %in% groups & sample_site == 'CE') |>
    group_by(sample_site, date) |> 
    summarize(
      pgC_L = sum(pgC_L)
    ) |> 
    mutate(month = month(date)) |> 
    group_by(sample_site, month) |> 
    summarize(
      mean_C = mean(pgC_L,na.rm = T),
      sd_C = sd(pgC_L, na.rm = T)
    )
  p = ggplot(data) +
    geom_point(
      aes(
        x = month,
        y = mean_C,
        color = sample_site
      )
    ) +
    geom_ribbon(
      aes(
        x = month,
        ymax = mean_C + sd_C,
        ymin = ifelse(
          (mean_C - sd_C) > 0,
          (mean_C - sd_C),
          0
        ),
        fill = sample_site
      ),
      alpha = 0.5
    )+
    theme_minimal()
  return(p)
} 

seasonality_plot(etx$names$diatom)
seasonality_plot(etx$names$mz)

# region \- environmental plots ------------------

etx$conc |> 
  filter(taxa %in% etx$names$diatom) |> 
  group_by(sample_site, date) |> 
  summarize(
    pgC_L = sum(pgC_L),
    sal = mean(sal)
  ) |> 
ggplot() +
  geom_point(
    aes(x = sal, y = pgC_L)
  ) +
  theme_minimal()

etx$conc |> 
  filter(taxa %in% c('Choreotrichia', 'Choreotrichia X'), 
         sample_site %in% c("CE", "CW"),
         month(yearmo) %in% c(06,07,08,09)) |> 
  group_by(sample_site, date) |> 
  summarize(
    pgC_L = sum(pgC_L),
    sal = mean(sal)
  ) |> 
ggplot() +
  geom_point(
    aes(x = sal, y = pgC_L)
  ) +
  theme_minimal()
