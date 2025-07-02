rm(list = ls())
library(ggplot2)
library(lubridate)
library(dplyr)
library(ggpubr)

# region Get Date Range -----------------------

etx <- readRDS('./data/02-full_merged.rds')

date_seq <- seq.Date(min(etx$conc$yearmo), max(etx$conc$yearmo), by = '1 month')

rm(etx)
# endregion -----------------------



# region \- environmental data --------------------------


environ <- readRDS("./data/01-environ_clean.rds")

environ$wind_avg$yearmo <- make_yearmo(environ$wind_avg)
environ$wind_avg <- environ$wind_avg |> filter(yearmo %in% date_seq)

wind <- environ$wind_avg |> 
  group_by(yearmo) |> 
  summarize(
    wind = mean(windspeed, na.rm = TRUE),
    PAR = mean(TotalPAR, na.rm = TRUE)
  )


# water quality -
environ$wq_sum$yearmo <- make_yearmo(environ$wq_sum)
environ$wq_sum <- environ$wq_sum |> filter(environ$wq_sum$yearmo %in% date_seq)

wq_sum <- environ$wq_sum |> 
  group_by(yearmo, sample_site) |> 
  summarize(
    temp = mean(t, na.rm = TRUE),
    sal = mean(sal, na.rm = TRUE),
  )


# chla lab-based measurement
environ$nut_avg$yearmo <- as.Date(
  paste(year(environ$nut_avg$date), month(environ$nut_avg$date), '01', sep = '-')
)

environ$nut_avg <- environ$nut_avg |> 
  filter(yearmo %in% date_seq)

chla_avg <- environ$nut_avg |> 
  group_by(yearmo, sample_site) |> 
  summarize(chl  = mean(CHLA_N, na.rm = TRUE))


# endregion

# region Plot TS -----------------------

ggplot(wind) +
  geom_line(
    aes(
      x = yearmo,
      y = wind
    )
  ) +
  theme_pubclean()

ggplot(chla_avg) +
  geom_line(
    aes(
      x = yearmo,
      y = chl,
      color = sample_site
    )
  ) +
  theme_pubclean()

ggplot() +
  geom_line(
    data = wq_sum,
    aes(
      x = yearmo,
      y = sal,
      color = sample_site
    )
  ) +
  # geom_rug(
  #   aes(x = discharge$Date, color = regime),
  #   sides = "b",
  # ) +
  theme_pubclean()

ggplot(wq_sum) +
  geom_line(
    aes(
      x = yearmo,
      y = temp,
      color = sample_site
    )
  ) +
  theme_pubclean()

ggplot(wq_sum) +
  geom_line(
    aes(
      x = yearmo,
      y = do,
      color = sample_site
    )
  ) +
  theme_pubclean()

ggplot(wq_sum) +
  geom_line(
    aes(
      x = yearmo,
      y = turb,
      color = sample_site
    )
  ) +
  theme_pubclean()



# endregion -----------------------




# region Supplemental Seasonality -----------------------

# possibly better to do averages and look at it all that way
plot(
  y = environ$wq_sum$t,
  x = yday(environ$wq$date),
  col = as.factor(environ$wq_sum$sample_site),
  pch = 16
)

plot(
  y = environ$wq_sum$DO,
  x = yday(environ$wq$date),
  col = as.factor(environ$wq_sum$sample_site),
  pch = 16
)

plot(
  y = environ$wq_sum$sal,
  x = yday(environ$wq$date),
  col = as.factor(environ$wq_sum$sample_site),
  pch = 16
)


plot(
  x = yday(environ$nut_avg$date),
  y = environ$nut_avg$CHLA_N,
  col = as.factor(environ$nut_avg$sample_site),
  pch = 16
)


plot(
  x = yday(environ$wind_avg$date),
  y = environ$wind_avg$windspeed,
  pch = 16
)

plot(
  x = yday(environ$wind_avg$date),
  y = environ$wind_avg$TotalPAR,
  pch = 16
)


# endregion -----------------------
