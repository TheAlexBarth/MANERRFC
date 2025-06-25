####
# Import Flow Gauge Status ######
####
rm(list = ls())
library(ggplot2)
library(lubridate)
library(dplyr)
library(tidyr)

# data are available from USGS river gauge
# publically available.
box_path <- '~/Library/CloudStorage/Box-Box/TGCRC Plankton Food Webs/Data'

ar <- paste0(
  box_path, '/USGS-GAUGE-AransasRv.csv'
) |> read.csv() |> 
  mutate(river = 'Aransas')

mr <- paste0(
  box_path, '/USGS-GAUGE-MissionRv.csv'
) |> read.csv() |> 
  mutate(river = 'Mission')

flow <- rbind(ar,mr)

flow$Date <- flow$Date |> as.Date()
flow$yearmo <- paste0(year(flow$Date), '-',month(flow$Date), '-01') |> 
  as.Date()

discharge <- flow |> 
  filter(year(Date) %in% c(2014:2021))

discharge <- discharge |> 
  select(-Site) |> 
  pivot_wider(values_from = Discharge, names_from = river)

cor(discharge$Mission, discharge$Aransas)

# pca for riverflow

river_pca <- prcomp(discharge[,c("Mission",'Aransas')])
summary(river_pca)


# define rolling window -------
# all(seq.Date(as.Date('2014-01-01'),as.Date('2021-12-31'),'1 day') %in% river_mat$Date)

high_flow <- river_pca$x[,1] > quantile(river_pca$x[,1], 0.95)


roll_window <- function(vect, k) {
  n <- length(vect)
  result <- rep(NA, n)
  half_k <- floor(k / 2)

  for (i in seq_len(n)) {
    start <- max(1, i - half_k)
    end <- min(n, i + half_k)
    result[i] <- any(vect[start:end], na.rm = TRUE)
  }

  return(result)
}

flow_cluster <- roll_window(high_flow, 120)

regime <- rep(NA, length(flow_cluster))
regime[flow_cluster] <- 'Wet'
regime[!flow_cluster] <- 'Dry'

ggplot() +
  geom_line(
    data = discharge,
    aes(
      x = as.Date(Date),
      y = Mission,
    ),
    color = 'blue',
    alpha = 0.75,size = .5
  ) +
  geom_line(
    data = discharge,
    aes(
      x = as.Date(Date),
      y = Aransas,
    ),
    color = 'red',
    alpha = 0.75,size = .5
  ) +  
  geom_rug(
    aes(x = discharge$Date, color = regime),
    sides = "b",
  ) +
  theme_minimal()
