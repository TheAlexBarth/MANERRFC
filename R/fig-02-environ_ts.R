rm(list = ls())
library(ggplot2)
library(lubridate)
library(dplyr)
library(ggpubr)
library(cowplot)
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
si <- readRDS('./data/01e-silicate.rds')
si <- si[which(year(si$yearmo) %in% c(2014:2021)),]
wind <- readRDS('./data/01c-wind_score.RDS')
wind <- wind[which(year(wind$yearmo) %in% c(2014:2021)),]
regime <- readRDS('./data/01d-river_regime.RDS')
regime <- regime[which(year(regime$yearmo) %in% c(2014:2021)),]

# join SiOH onto the continuous monthly nutrient record (not gated to
# phytoplankton sampling occasions) and add total N, molar (uM)
# concentrations, and N:P / N:Si ratios (see compute_molar_ratios() in utils.R)
nut <- swmp |>
  left_join(si, by = c('sampling_site', 'yearmo'))
nut <- compute_molar_ratios(nut)

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

swmp_plotter <- function(var, ylab, legend = FALSE) {
  p = ggplot(nut) +
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
  if (!legend) p <- p + theme(legend.position = 'none')
  return(p)
}

temp_plot <- swmp_plotter('temp', 'Temperature [deg.C]')
sal_plot <- swmp_plotter('sal', 'Salinity')
P_plot <- swmp_plotter('P', 'PO4 [mg/L]')
NH4_plot <- swmp_plotter('NH4', 'NH4 [mg/L]')
N_plot <- swmp_plotter('N', 'NO23 [mg/L]')
Si_plot <- swmp_plotter('SiOH', 'SiOH [mg/L]')
NP_plot <- swmp_plotter('NP_ratio', 'N:P [uM/uM]')
NSi_plot <- swmp_plotter('NSi_ratio', 'N:Si [uM/uM]')

# endregion


# region \- full plot -------------------------------------
# regime + wind each keep their own legend (different series, not site);
# the site-colored panels (temp through N:Si) share one common site legend,
# extracted once and placed alongside the stacked panels, rather than
# repeating it 8 times.
#
# NOTE: nesting two ggarrange()/plot_grid() calls (e.g. an outer arrangement
# of two already-arranged blocks, one carrying a ggarrange(common.legend=
# TRUE) legend) corrupts the render - panels come out as solid black boxes.
# Extracting the legend once and building a single, flat plot_grid() avoids
# the nested-arrangement grobs that trigger it.

margin_theme <- theme(plot.margin = unit(c(0, .2, 0, 0.2), 'lines'))

site_legend <- get_legend(
  swmp_plotter('temp', 'Temperature [deg.C]', legend = TRUE) +
    theme(legend.position = 'right')
)

main_grid <- plot_grid(
  regime_plot + margin_theme,
  wind_plot + margin_theme,
  temp_plot + margin_theme,
  sal_plot + margin_theme,
  P_plot + margin_theme,
  NH4_plot + margin_theme,
  N_plot + margin_theme,
  Si_plot + margin_theme,
  NP_plot + margin_theme,
  NSi_plot + margin_theme,
  ncol = 1,
  align = 'v',
  labels = LETTERS[1:10]
)

full_plot <- plot_grid(main_grid, site_legend, ncol = 2, rel_widths = c(10, 1.4))

ggsave(
  './output/fig02-environ.pdf', full_plot,
  height = 340, width = 190,
  units = 'mm', dpi = 600, bg = 'white'
)

# endregion

# region \- ratio trend table -----------------------
# per-site lm(ratio ~ yearmo) for the N:P / N:Si panels above: slope is
# change in ratio per day, multiplied by 365.25 for an annual rate.
trend_summary <- function(var) {
  nut |>
    filter(!is.na(.data[[var]])) |>
    group_by(sampling_site) |>
    summarize(
      n = n(),
      slope_per_yr = tryCatch(
        coef(lm(.data[[var]] ~ as.numeric(yearmo)))[2] * 365.25,
        error = function(e) NA_real_
      ),
      p_value = tryCatch(
        summary(lm(.data[[var]] ~ as.numeric(yearmo)))$coefficients[2, 4],
        error = function(e) NA_real_
      ),
      .groups = 'drop'
    ) |>
    mutate(ratio = var)
}

trend_table <- bind_rows(
  trend_summary('NP_ratio'),
  trend_summary('NSi_ratio')
)

write.csv(trend_table, './output/fig02-nutrient_ratio_trends.csv', row.names = FALSE)
# endregion -----------------------
