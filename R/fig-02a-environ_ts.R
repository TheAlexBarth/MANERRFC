######
# Figure 2a - Environmental time series
#
# Physical / hydrographic drivers: river inflow (Aransas + Mission discharge
# with wet/dry regime), wind forcing, water temperature, and salinity.
# Companion to fig-02b (nutrients). Single-column width (89 mm ~ 3.5"),
# panels stacked vertically with a shared x-axis.
######

rm(list = ls())
library(ggplot2)
library(lubridate)
library(dplyr)
library(ggpubr)
library(cowplot)
source('./R/utils.R')

# region \- data prep --------------------------
yrs <- 2014:2021

swmp <- readRDS('./data/01a-swmp_wq_data.rds')
swmp <- swmp[year(swmp$yearmo) %in% yrs, ]
swmp$sampling_site <- factor(swmp$sampling_site, levels(site_factors))

wind <- readRDS('./data/01c-wind_score.RDS')
wind <- wind[year(wind$yearmo) %in% yrs, ]

regime <- readRDS('./data/01d-river_regime.RDS')
regime <- regime[year(regime$yearmo) %in% yrs, ]

# shared x-axis so every stacked panel lines up on the same monthly grid;
# 2-year breaks keep the labels legible at single-column width
x_axis <- scale_x_date(
  limits = as.Date(c('2013-12-15', '2022-01-15')),
  date_breaks = '2 years', date_labels = '%Y', expand = c(0, 0)
)

base_theme <- theme_pubclean(base_size = 8) +
  theme(
    plot.margin = margin(1, 5, 1, 4),
    axis.title.y = element_text(size = 7.5),
    axis.text = element_text(size = 6.5),
    legend.text = element_text(size = 6),
    legend.key.size = unit(8, 'pt'),
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(2, 'pt')
  )

# drop the x tick labels on the upper panels (bottom panel carries them)
strip_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# endregion

# region \- panels --------------------------
# discharge is plotted in thousands of ft^3/s so its y tick labels are ~2 digits
# like the other panels - this keeps the y-axis titles aligned down the column
# (the wide 5-digit labels were pushing this title out of line)
inflow_depth <- 0.06 * max(c(regime$Aransas, regime$Mission), na.rm = TRUE) / 1000
inflow_runs <- transform(regime_run_bands(regime), ymin = -inflow_depth, ymax = 0)

inflow_plot <- ggplot(regime) +
  geom_line(aes(Date, Aransas / 1000, color = 'Aransas')) +
  geom_line(aes(Date, Mission / 1000, color = 'Mission')) +
  geom_rect(
    data = inflow_runs,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = regime),
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = c(Aransas = '#01086b', Mission = '#979df7')) +
  scale_fill_manual(values = regime_cols) +
  labs(x = '', y = expression('Discharge [10'^3~'ft'^3~'s'^-1 * ']'),
    color = '', fill = '') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(-inflow_depth, NA)) +
  x_axis + base_theme

wind_plot <- ggplot(wind, aes(x = yearmo)) +
  geom_line(aes(y = ptat_wind_ms, linetype = 'PTAT')) +
  geom_line(aes(y = rcpt_wind_ms, linetype = 'RCPT')) +
  geom_line(aes(y = awrt_wind_ms, linetype = 'AWRT')) +
  geom_line(aes(y = wind_pca, linetype = 'PC1')) +
  scale_linetype_manual(values = c(
    PTAT = 'twodash', RCPT = 'dashed', AWRT = 'dotdash', PC1 = 'solid'
  )) +
  labs(x = '', y = expression('Wind speed [m s'^-1 * ']'), linetype = '') +
  x_axis + base_theme

site_ts <- function(var, ylab) {
  ggplot(swmp) +
    geom_line(aes(x = yearmo, y = .data[[var]], color = sampling_site)) +
    scale_color_manual(values = site_cols, drop = FALSE) +
    labs(x = '', y = ylab, color = NULL) +
    x_axis + base_theme + theme(legend.position = 'none')
}

temp_plot <- site_ts('temp', expression('Temperature [' * degree * 'C]'))
sal_plot  <- site_ts('sal', 'Salinity')
# endregion

# region \- assemble --------------------------
# Legends live in a right-hand column, each aligned to the panel(s) it
# describes: inflow and wind get their own; temp + sal share one title-less
# site legend centered across both. Panels themselves are legend-free so they
# keep identical widths and their x-axes stay aligned down the column.
inflow_leg <- get_legend(
  inflow_plot + theme(legend.position = 'right') +
    guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1))
)
wind_leg <- get_legend(
  wind_plot + theme(legend.position = 'right') +
    guides(linetype = guide_legend(ncol = 1))
)
site_leg <- get_legend(
  site_ts('temp', 'x') + theme(legend.position = 'right') +
    guides(color = guide_legend(ncol = 1))
)

# theme_pubclean() defaults to a top legend, so explicitly drop it on the
# inflow/wind panels (temp/sal already set legend.position = 'none')
no_leg <- theme(legend.position = 'none')

panels <- plot_grid(
  inflow_plot + strip_x + no_leg,
  wind_plot + strip_x + no_leg,
  temp_plot + strip_x,
  sal_plot,
  ncol = 1, align = 'v', axis = 'lr',
  labels = c('A', 'B', 'C', 'D'), label_size = 9,
  rel_heights = c(1, 1, 1, 1.13)  # bottom panel carries the x-axis
)

# legend column: heights mirror the panel rows so each legend sits beside its
# panel; the site legend spans the temp + sal rows (1 + 1.13)
legend_col <- plot_grid(
  inflow_leg, wind_leg, site_leg,
  ncol = 1, rel_heights = c(1, 1, 2.13)
)

fig_2a <- plot_grid(panels, legend_col, ncol = 2, rel_widths = c(3.6, 1))

ggsave('./output/fig02a-environ.pdf', fig_2a,
  width = 89, height = 190, units = 'mm', dpi = 600, bg = 'white')
# endregion
