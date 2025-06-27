rm(list = ls())
library(dplyr)
library(ggplot2)
library(ggpubr)

years <- 2011:2021

# region PTAT2 -----------------------
# https://www.ndbc.noaa.gov/station_history.php?station=ptat2

ptat2 <- list()

for(year in years) {
  url = paste0("https://www.ndbc.noaa.gov/view_text_file.php?filename=ptat2c",year,".txt.gz&dir=data/historical/cwind/")
  raw = read.table(url)
  names(raw) <- c("Year","month","day","hour",'minute', "wind_dir",'windspeed_ms','GDR','GST','GTIME')
  raw$windspeed_ms[raw$windspeed_ms >= 99.0] <- NA
  raw$yearmo <- paste(
    raw$Year, raw$month, '01',sep = '-'
  ) |> as.Date()
  ptat2[[year]] <- raw |> 
    select(
      yearmo,
      windspeed_ms
    ) |> 
    group_by(
      yearmo
    ) |> 
    summarize(ptat_wind_ms = mean(windspeed_ms,na.rm = TRUE))
}

ptat2 <- do.call(rbind, ptat2)

# endregion -----------------------

# region Section -----------------------
#https://www.ndbc.noaa.gov/view_text_file.php?filename=rcpt2h2014.txt.gz&dir=data/historical/stdmet/

rcpt2 <- list()
for(year in years) {
  url = paste0("https://www.ndbc.noaa.gov/view_text_file.php?filename=rcpt2h",year,".txt.gz&dir=data/historical/stdmet/")
  raw = read.table(url)
  names(raw)[1:7] <- c("Year","month","day","hour",'minute', "wind_dir",'windspeed_ms')
  raw$windspeed_ms[raw$windspeed_ms >= 99.0] <- NA
  raw$yearmo <- paste(
    raw$Year, raw$month, '01',sep = '-'
  ) |> as.Date()
  rcpt2[[year]] <- raw |> 
    select(
      yearmo,
      windspeed_ms
    ) |> 
    group_by(
      yearmo
    ) |> 
    summarize(rcpt_wind_ms = mean(windspeed_ms,na.rm = TRUE))
}
rcpt2 <- do.call(rbind, rcpt2)

# endregion -----------------------

# region Section -----------------------

awrt2h <- list()
for(year in years) {
  url = paste0("https://www.ndbc.noaa.gov/view_text_file.php?filename=awrt2h",year,".txt.gz&dir=data/historical/stdmet/")
  raw = tryCatch(read.table(url), error = function(e) return(NULL))
  if(is.null(raw)) {next}
  names(raw)[1:7] <- c("Year","month","day","hour",'minute', "wind_dir",'windspeed_ms')
  raw$windspeed_ms[raw$windspeed_ms >= 99.0] <- NA
  raw$yearmo <- paste(
    raw$Year, raw$month, '01',sep = '-'
  ) |> as.Date()
  awrt2h[[year]] <- raw |> 
    select(
      yearmo,
      windspeed_ms
    ) |> 
    group_by(
      yearmo
    ) |> 
    summarize(awrt_wind_ms = mean(windspeed_ms,na.rm = TRUE))
}
awrt2h <- do.call(rbind, awrt2h)

# endregion -----------------------

full_wind <- ptat2 |> 
  full_join(
    rcpt2, by = 'yearmo'
  ) |> 
  full_join(
    awrt2h, by = 'yearmo'
  ) |> 
  unique()


# ggplot(full_wind, aes(x = yearmo)) +
#   geom_line(aes(y = ptat_wind_ms), color = 1) +
#   geom_line(aes(y = rcpt_wind_ms), color = 2) +
#   geom_line(aes(y = awrt_wind_ms), color = 3) +
#   theme_pubclean()

cor(full_wind[,-1], use = 'na.or.complete')

nb_comp <- missMDA::estim_ncpPCA(full_wind[,-1])

imputed_wind <- missMDA::imputePCA(full_wind[,-1])

plot(imputed_wind$completeObs[,1], type = 'l', ylim = c(-3,11))
lines(imputed_wind$completeObs[,2])
lines(imputed_wind$completeObs[,3])

pca <- prcomp(imputed_wind$completeObs)
summary(pca)
lines(pca$x[,1])

outdf <- cbind(full_wind, `wind_pca` = pca$x[,1])
saveRDS(outdf, './data/01c-wind_score.RDS')
