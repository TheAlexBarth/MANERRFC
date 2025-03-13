######
# Environmental Formatting
######


rm(list = ls())
library(dplyr)
library(lubridate)

nut <- readRDS('./data/swmp_nut.rds')
nut$sample_site <- nut$StationCode |> substr(4,5) |> toupper()
wq <- readRDS('./data/swmp_wq.rds')
wq$sample_site <- wq$StationCode |> substr(4,5) |> toupper()
PDSI <- read.csv('./data/PDSI_South-texas.csv')

wind <- read.csv('./data/MARCEMET_R.csv')
wind <- wind |> mutate(
  DateTimeStamp = mdy_hm(DateTimeStamp) |> format("%Y/%d/%m")
)

# averaging values

nut_avg <- nut |> 
  mutate(date = as.Date(DateTimeStamp)) |> 
  select(date, PO4F, NH4F, NO23F, CHLA_N, sample_site) |> 
  group_by(date, sample_site) |> 
  summarize(
    P = mean(PO4F, na.rm = T),
    NH4 = mean(NH4F, na.rm = T),
    N = mean(NO23F, na.rm = T),
    CHLA_N = mean(CHLA_N, na.rm = T)
  ) |> 
  mutate(sample_id = paste(sample_site, date, sep = "_")) |> 
  ungroup()

# average daily wind speed values 

wind_avg <- wind |>
  mutate(date = as.Date(DateTimeStamp)) |>
  select(date, WSpd, MaxWSpd, TotPAR) |>
  group_by(date) |>
  summarize(
    windspeed = mean(WSpd, na.rm = T),
    max_windspeed = mean(MaxWSpd, na.rm = T),
    TotalPAR = mean(TotPAR, na.rm = T)
  ) |>
  mutate(sample_id = paste(date)) |>
  ungroup()

# water quality daily averages

wq_sum <- wq |> 
  mutate(date = as.Date(DateTimeStamp)) |> 
  group_by(date, sample_site) |> 
  summarize(
          t = mean(Temp, na.rm = T),
          t_max = max(Temp, na.rm = T),
          t_min = min(Temp, na.rm = T),
          s_max = max(Sal, na.rm = T),
          sal = mean(Sal, na.rm = T),
          DO = mean(DO_mgl, na.rm = T),
          DO_min = min(DO_mgl, na.rm = T),
          Chl = mean(ChlFluor, na.rm = T),
          Chl_max = max(ChlFluor, na.rm = T),
          Chl_min = min(ChlFluor, na.rm = T)
          ) |> 
  mutate(sample_id = paste(sample_site, date, sep = '_')) |> 
  ungroup()

# region \- clean up data -------------------

extreme_to_na <- function(vect, sd_away) {
  vect[which(is.infinite(vect))] <- NA
  mv = mean(vect, na.rm = T)
  sv = sd(vect, na.rm = T)

  vout <- vect |> 
    sapply(
      function(x)
      ifelse(x > mv+sd_away*sv, NA, x)
    )
  
  return(vout)
}

# there's some crazy values out there so I'm just removing massive outliers
# that snuck through the QAQC

wq_sum <- wq_sum |> 
  mutate(
    across(c(s_max, sal, DO, DO_min, Chl, Chl_max, Chl_min), ~ extreme_to_na(.,5))
  )



saveRDS(
  list(
    nut_avg = nut_avg,
    wq_sum = wq_sum,
    PDSI = PDSI,
    wind_avg = wind_avg
  ),
  './data/01-environ_clean.rds'
)
