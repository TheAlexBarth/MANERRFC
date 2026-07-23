rm(list = ls())
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(lubridate)

source('./R/utils.R')


chl <- readRDS('./data/01b-size_frac_chl.RDS')
regime <- readRDS("./data/01d-river_regime.RDS")

# wet/dry regime drawn as a continuous band below y = 0 (contiguous run
# rectangles, see regime_run_bands() in utils.R) instead of a dense daily rug
ymax_chl <- 30
regime_depth <- 0.06 * ymax_chl
regime_runs <- transform(regime_run_bands(regime), ymin = -regime_depth, ymax = 0)

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
    geom_rect(data = subset(regime_runs, regime == 'Dry'),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = regime_cols[['Dry']], inherit.aes = FALSE) +
    geom_rect(data = subset(regime_runs, regime == 'Wet'),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = regime_cols[['Wet']], inherit.aes = FALSE) +
    scale_fill_manual(values = size_cols, labels = size_labels)+
    labs(
      x = "", y = expression(paste("Chl-a ", "[",mu * g~L^{-1},"]")),
      fill = "", subtitle = site_labels[[loc]]
    )+
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
    coord_cartesian(ylim = c(-regime_depth, ymax_chl))+
    theme_pubclean(base_size = 8)
  return(plt)
}

all_plots <- list()
for(site in names(site_cols)) {
  all_plots[[site]] <- site_chl_plot(site)
}

# panels stacked in one column, but the figure is sized to two-column print width
outplot <- ggarrange(
  plotlist = all_plots,
  ncol = 1,
  common.legend = TRUE,
  legend = 'bottom',
  labels = LETTERS
)

ggsave(
  './output/fig03-chla_plot.pdf',outplot,
  width = 178, height = 200, units = "mm",
  dpi = 600
)
