rm(list = ls())
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

source('./R/utils.R')


chl <- readRDS('./data/01b-size_frac_chl.RDS')
regime <- readRDS("./data/01d-river_regime.RDS")


site_chl_plot <- function(loc) {
  sub_data <- chl |> 
    filter(site == loc) |> 
    pivot_longer(
      cols = c('micro','nano','pico'),
      names_to = 'frac',
      values_to = 'ugChl_L'
    )
  
  plt <- ggplot() +
    geom_area(
      data = sub_data,
      aes(
        x = yearmo,
        y = ugChl_L,
        fill = frac
      ),
      position = 'stack'
    ) + 
    geom_rug(
      data = regime |> 
        filter(year(Date) %in% c(min(sub_data$year):max(sub_data$year))),
      aes(
        x = Date,
        color = regime
      )
    )+
    scale_color_manual(values = regime_cols)+
    scale_fill_manual(values = size_cols)+
    guides(color = 'none')+
    labs(
      x = "", y = expression(paste("Chl-a ", "[",mu * g~L^{-1},"]")), 
      fill = "", subtitle = loc
    )+
    theme_pubclean(base_size = 8)
  return(plt)
}

all_plots <- list()
for(site in names(site_cols)) {
  all_plots[[site]] <- site_chl_plot(site)
}

outplot <- ggarrange(
  plotlist = all_plots,
  ncol = 1,
  common.legend = TRUE,
  legend = 'bottom'
)

ggsave(
  './output/fig03-chla_plot.pdf',outplot,
  width = 175, height = 175, units = "mm",
  dpi = 600
)
