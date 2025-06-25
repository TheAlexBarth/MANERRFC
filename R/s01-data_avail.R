rm(list = ls())
library(lubridate)

etx <- readRDS('./data/02-full_merged.rds')


# data availability plot

png('./output/s01-data_avail.png')
plot(
  x = month(etx$conc$yearmo),
  y = year(etx$con$yearmo),
  xlab = 'Month', ylab = 'Year'
)
dev.off()
