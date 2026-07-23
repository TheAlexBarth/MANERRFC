######
# Figure 2b - Nutrient time series
#
# Monthly SWMP nutrient concentrations (PO4, NH4, NO23, SiOH) and molar
# stoichiometric ratios (N:P, N:Si) by site. Companion to fig-02a
# (environmental drivers). Full manuscript width (178 mm ~ 7").
#
# Also writes the per-site N:P / N:Si linear trend table used in the text.
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

# add totalN, molar concentrations, and N:P / N:Si ratios (see utils.R)
nut <- compute_molar_ratios(swmp)

# wet/dry river regime, drawn as a continuous band BELOW y = 0 on every panel
# (matches fig-02a / fig-07): contiguous run-rectangles rather than a dense
# daily rug, in fixed colors so it doesn't collide with the site color scale on
# the lines. Its legend is built separately below.
regime <- readRDS('./data/01d-river_regime.RDS')
regime <- regime[year(regime$yearmo) %in% yrs, ]
regime_runs <- regime_run_bands(regime)

# shared x-axis so every stacked panel lines up on the same monthly grid;
# odd-year labels only, so the left column's last tick can't crowd the right
x_axis <- scale_x_date(
  limits = as.Date(c('2013-12-15', '2022-01-15')),
  breaks = seq(as.Date('2015-01-01'), as.Date('2021-01-01'), by = '2 years'),
  date_labels = '%Y', expand = c(0, 0)
)

base_theme <- theme_pubclean(base_size = 8) +
  theme(
    plot.margin = margin(2, 10, 2, 4),
    axis.title.y = element_text(size = 7.5),
    axis.text = element_text(size = 6.5),
    legend.text = element_text(size = 6.5),
    legend.key.size = unit(9, 'pt'),
    legend.margin = margin(0, 0, 0, 0)
  )
# endregion

# region \- panels --------------------------
nutrient_ts <- function(var, ylab) {
  depth <- 0.06 * max(nut[[var]], na.rm = TRUE)
  if (!is.finite(depth) || depth <= 0) depth <- 1
  runs <- transform(regime_runs, ymin = -depth, ymax = 0)
  ggplot() +
    geom_line(data = nut, aes(x = yearmo, y = .data[[var]], color = sampling_site)) +
    geom_rect(data = subset(runs, regime == 'Dry'),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = regime_cols[['Dry']], inherit.aes = FALSE) +
    geom_rect(data = subset(runs, regime == 'Wet'),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = regime_cols[['Wet']], inherit.aes = FALSE) +
    scale_color_manual(values = site_cols, drop = FALSE) +
    labs(x = '', y = ylab, color = "") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(ylim = c(-depth, NA)) +
    x_axis + base_theme + theme(legend.position = 'none')
}

P_plot   <- nutrient_ts('P',    expression(PO[4] ~ '[mg L'^-1 * ']'))
NH4_plot <- nutrient_ts('NH4',  expression(NH[4] ~ '[mg L'^-1 * ']'))
N_plot   <- nutrient_ts('N',    expression(NO[2] * '+' * NO[3] ~ '[mg L'^-1 * ']'))
Si_plot  <- nutrient_ts('SiOH', expression(SiOH ~ '[mg L'^-1 * ']'))
NP_plot  <- nutrient_ts('NP_ratio',  'N:P')
NSi_plot <- nutrient_ts('NSi_ratio', 'N:Si')
# endregion

# region \- assemble --------------------------
# all six panels are site-colored (one shared site legend) plus the wet/dry rug
# (its own legend, built from a dummy fill scale). Both sit at the bottom.
# Single flat plot_grid(grid, legend) - see fig-02a note on nested legends.
site_legend <- get_legend(
  nutrient_ts('P', 'x') +
    theme(legend.position = 'bottom') +
    guides(color = guide_legend(nrow = 1))
)

regime_legend <- get_legend(
  ggplot(regime, aes(x = Date, y = 1, fill = regime)) +
    geom_tile() +
    scale_fill_manual('', values = regime_cols) +
    base_theme +
    theme(legend.position = 'bottom') +
    guides(fill = guide_legend(nrow = 1))
)

legends <- plot_grid(site_legend, regime_legend, nrow = 1, rel_widths = c(1, 0.7))

grid_2b <- plot_grid(
  P_plot, NH4_plot,
  N_plot, Si_plot,
  NP_plot, NSi_plot,
  ncol = 2, align = 'hv', axis = 'tblr',
  labels = 'AUTO', label_size = 9
)

fig_2b <- plot_grid(grid_2b, legends, ncol = 1, rel_heights = c(1, 0.05))

ggsave('./output/fig02b-nutrients.pdf', fig_2b,
  width = 178, height = 160, units = 'mm', dpi = 600, bg = 'white')
# endregion

# region \- ratio trend table --------------------------
# per-site lm(ratio ~ yearmo): slope is change per day, x365.25 for annual rate
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
# endregion
