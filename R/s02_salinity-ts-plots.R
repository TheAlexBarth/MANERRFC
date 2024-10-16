####
# Quick Salinity Plots #####
####


rm(list = ls()) # this wipes your environment
library(EcotaxaTools) # this is my package
library(lubridate)
library(dplyr)
library(ggplot2)
library(tidyr)


wq <- readRDS('./data/01_swmp_wq.rds')
wq$sample_site <- wq$StationCode |> substr(4,5) |> toupper()

wq <- wq |> 
  filter(sample_site %in% c('AB','CE','CW')) |> 
  filter(DateTimeStamp < as.Date('2023-01-01'))

wq$Date <- as.Date(wq$DateTimeStamp)

wq_daysum <- wq |> 
  group_by(Date, sample_site) |> 
  summarize(Temp = mean(Temp, na.rm = T),
            Sal = mean(Sal, na.rm = T),
            Depth = mean(Depth, na.rm = T),
            pH = mean(pH, na.rm = T),
            ChlFluor = mean(ChlFluor, na.rm = T))


###
# Time Series Plot #####
###

ggplot(wq_daysum) +
  geom_line(
    aes(x = Date, y = Sal, color = sample_site)
  ) + 
  labs(x = "", y = "Salinity") +
  scale_color_manual(values = gg_cbb_col(3)) + 
  theme_minimal() +
  theme()


####
# Mo-day ####
###
wq <- wq |> 
  mutate(mo_day = as.Date(
    paste('2000',month(Date), day(Date),sep = '-')
  ))

wq_daysum$mo_day = as.Date(
  paste('2000',month(wq_daysum$Date), day(wq_daysum$Date),sep = '-')
)

wq_moday <- wq |> 
  group_by(mo_day, sample_site) |> 
  summarize(
    sal = mean(Sal, na.rm = T),
    sd_sal = sd(Sal, na.rm = T)
  )


point_plot <- ggplot() +
  geom_point(
    data = wq_daysum,
    aes(x = mo_day, y = Sal, color = sample_site),
    size = 0.25, alpha = 0.1
  )+
  geom_smooth(
    data = wq_moday,
    aes(x = mo_day, y = sal, color = sample_site),
    size = 1.5, span = 0.1
  )+
  scale_color_manual(values = gg_cbb_col(3)) +
  labs(x = "", y = "Salinity")+
  scale_x_date(date_labels = '%b')+
  theme_minimal() +
  theme(legend.position = 'none', axis.title = element_text(size = 8),
axis.text = element_text(size = 6))

ggsave('./output/s02b_salinity-points.pdf', point_plot,
height = 3, width = 6, units = 'in', dpi = 600)
