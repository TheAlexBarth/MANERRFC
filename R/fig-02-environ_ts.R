rm(list = ls())
library(ggplot2)
library(lubridate)
library(dplyr)
library(ggpubr)
source('./R/utils.R')

# region Get Date Range -----------------------

etx <- readRDS('./data/02-full_merged.rds')

date_seq <- seq.Date(min(etx$conc$yearmo), max(etx$conc$yearmo), by = '1 month')

rm(etx)
# endregion -----------------------



# region \- prep data data --------------------------

swmp <- readRDS("./data/01a-swmp_wq_data.rds")
swmp <- swmp[which(year(swmp$yearmo) %in% c(2014:2021)),]
swmp$sampling_site <- factor(swmp$sampling_site, levels(site_factors))
wind <- readRDS('./data/01c-wind_score.RDS')
wind <- wind[which(year(wind$yearmo) %in% c(2014:2021)),]
regime <- readRDS('./data/01d-river_regime.RDS')
regime <- regime[which(year(regime$yearmo) %in% c(2014:2021)),]


# endregion

# region Plot TS -----------------------

regime_plot <- ggplot(regime) +
  geom_line(
    aes(
      x = Date,
      y = Aransas,
      color = "Aransas"
    ),
  ) +
  geom_line(
    aes(
      x = Date,
      y = Mission,
      color = 'Mission'
    ),
  ) +
  geom_rug(
    aes(
      x = Date,
      color = regime
    )
  )+
  scale_color_manual(values = c(
    `Aransas` = "#01086b", `Mission` = '#979df7',
    regime_cols
  )) +
  labs(x = "", y = "Discharge [CuFt/s]", color = "")+
  theme_pubclean(base_size = 7)+
  theme(legend.position = 'top')

# endregion -----------------------


# region \- wind---------------

wind_plot <- ggplot() +
  geom_line(
    data = wind,
    aes(
      x = yearmo,
      y = ptat_wind_ms,
      linetype = 'PTAT'
    )
  ) +
  geom_line(
    data = wind,
    aes(
      x = yearmo,
      y = rcpt_wind_ms,
      linetype = 'RCPT'
    )
  ) +  
  geom_line(
    data = wind,
    aes(
      x = yearmo,
      y = awrt_wind_ms,
      linetype = 'AWRT'
    )
  ) +
  geom_line(
    data = wind,
    aes(
      x = yearmo,
      y = wind_pca,
      linetype = 'Wind Score (1st Eigenvalue)'
    )
  ) +
  scale_linetype_manual(
    values = c(
      `PTAT` = 'twodash',`RCPT` = 'dashed', `AWRT` = "dotdash",
      `Wind Score (1st Eigenvalue)` = 'solid'
    )
  ) +
  labs(x = "", y = 'Wind Speed [m/s]', linetype = "")+ 
  theme_pubclean(base_size = 7)
  
# endregion

# region \- swmp plots -----------------

swmp_plotter <- function(var, ylab) {
  p = ggplot(swmp) +
    geom_line(
      aes(
        x = yearmo,
        y = .data[[var]],
        color = sampling_site
      )
    )+
    scale_color_manual(values = site_cols) +
    labs(x = "", y = ylab, color = "")+
    theme_pubclean(base_size = 7)
  return(p)
}

temp_plot <- swmp_plotter('temp','Temperature [deg.C]')
sal_plot <- swmp_plotter('sal','Salinity')
P_plot <- swmp_plotter('P','PO4') 
NH4 <- swmp_plotter('NH4','NH4')
N <- swmp_plotter('N','NO23')

# endregion


# region \- full plot -------------------------------------

full_plot <- ggarrange(
  regime_plot + theme(plot.margin = unit(c(0,.2,0,0.2), 'lines')),
  wind_plot + theme(plot.margin = unit(c(0,.2,0,0.2), 'lines')),
  temp_plot + theme(plot.margin = unit(c(0,.2,0,0.2), 'lines')),
  sal_plot + theme(plot.margin = unit(c(0,.2,0,0.2), 'lines')),
  P_plot + theme(plot.margin = unit(c(0,.2,0,0.2), 'lines')),
  NH4 + theme(plot.margin = unit(c(0,.2,0,0.2), 'lines')),
  N + theme(plot.margin = unit(c(0,.2,0,0.2), 'lines')),
  ncol = 1,
  align = 'v'
)


ggsave(
  './output/fig02-environ.pdf',full_plot,
  height = 270, width = 170,
  units = 'mm', dpi = 600
)

# endregion
