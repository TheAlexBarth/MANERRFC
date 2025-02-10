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
  filter(sample_site %in% c('AB','CE','CW', 'SC')) |> 
  filter(DateTimeStamp < as.Date('2023-01-01'))

wq$Date <- as.Date(wq$DateTimeStamp)

wq_daysum <- wq |> 
  group_by(Date, sample_site) |> 
  summarize(Temp = mean(Temp, na.rm = T),
            Sal = mean(Sal, na.rm = T),
            Depth = mean(Depth, na.rm = T),
            pH = mean(pH, na.rm = T),
            ChlFluor = mean(ChlFluor, na.rm = T))


saveRDS(wq_daysum, './output/s02-temp-wq.rds')

###
# Time Series Plot #####
###

ggplot(wq_daysum |> 
    filter(sample_site %in% c('CE', 'CW')) |> 
    group_by(Date) |> 
    summarise(sal = mean(Sal))
    ) +
  geom_line(
    aes(x = Date, y = sal),
    linewidth = 3, color = '#dbf6f6'
  ) + 
  labs(x = "", y = "") +
  scale_color_manual(values = gg_cbb_col(3)) + 
  theme_minimal() +
  theme(
    axis.text = element_text(size = 28, color = 'white'),
    panel.grid = element_line(color = 'grey50'),
    plot.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.background = element_rect(fill = 'transparent')
  )
ggsave('./output/s02_sal-ts.pdf',
height = 6, width = 22, units = 'in')
  

# region \- chlorophyll reg ------

# ggplot(
#   wq_daysum |> 
#     filter(sample_site %in% c('CE', 'CW')) |> 
#     group_by(Date) |> 
#     summarise(chl = mean(ChlFluor), sal = mean(Sal))
#   ) +
#   geom_point(
#     aes(
#       x = sal,
#       y = chl
#     )
#   ) +
#   theme_minimal()


# chl_mod <- lm(
#   chl ~ sal,
#   data =  wq_daysum |> 
#     filter(sample_site %in% c('CE', 'CW')) |> 
#     group_by(Date) |> 
#     summarise(chl = mean(ChlFluor), sal = mean(Sal))
# )



# region \- chl ts 


ggplot(wq_daysum |> 
  filter(sample_site %in% c('CE', 'CW')) |> 
  group_by(Date) |> 
  summarise(chl = mean(ChlFluor))
  ) +
geom_line(
  aes(x = Date, y = chl),
  linewidth = 3, color = '#dbf6f6'
) + 
labs(x = "", y = "") +
scale_color_manual(values = gg_cbb_col(3)) + 
theme_minimal() +
theme(
  axis.text = element_text(size = 28, color = 'white'),
  panel.grid = element_line(color = 'grey50'),
  plot.background = element_rect(fill = 'transparent', color = 'transparent'),
  panel.background = element_rect(fill = 'transparent')
)


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

# MARK:

ggplot(data = wq_moday |> 
  filter(sample_site == 'SC' &
    mo_day > as.Date("2000-06-01") & mo_day < as.Date("2000-11-01")
  )) +
  geom_line(
    aes(
      x = mo_day,
      y = sal
    )
  ) +
  geom_ribbon(
    aes(
      x = mo_day,
      ymin = sal-sd_sal,
      ymax = sal+sd_sal
    ),
    color = '#56565600', alpha = 0.25
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
