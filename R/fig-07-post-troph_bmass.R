rm(list = ls())
library(ggplot2)
library(ggpubr)
library(dplyr)
library(lubridate)
source('./R/utils.R')

# NOTE: the "non-chlorophyll" predictor block for autotrophs (intercept,
# wind, temp, sal, P, NH4, N, SiOH - i.e. everything except chl) is sized
# dynamically via `1:(ncol(post$preds) - 3)` rather than hardcoded, since
# stan/count_mod.stan always puts chl (micro/nano/pico) in the final 3
# columns regardless of how many nutrient predictors precede them.
#
# NOTE: biomass is plotted via `median` rather than `mean`. mean_bmass is
# back-transformed from a log-normal (exp(eta_r + sigma_r^2/2)), which is
# extremely sensitive to the tail of sigma_r - a handful of MCMC draws with
# large sigma_r can blow the arithmetic mean up by orders of magnitude while
# barely nudging the 95% quantile ribbon. The median is robust to that and
# stays consistent with the ribbon (matches fig-06's `mid = median(est)`).
post <- readRDS("./data/04-post_micro.RDS")

conc_raw <- readRDS("./data/02-full_merged.rds")$conc
regime <- readRDS('./data/01d-river_regime.rds')
post$alpha |> dim()

nonchl_cols <- 1:(ncol(post$preds) - 3)

# region CALCULATE Lambda exp


troph_list <- list()
for(troph in troph_factors) {
  troph_list[[troph]] <- list()
  troph_idx <- which(levels(troph_factors) == troph)

  mean_bmass <- exp(post$eta_r[,troph_idx] + (post$sigma_r[,troph_idx]^2/2)) # pgC

  for(site in site_factors) {
    site_idx <- which(levels(site_factors) == site)

    data_idx <- which(conc_raw$functional_role == troph & conc_raw$sample_site == site)
    if(troph_idx < 3) {
      lambda <- as.matrix(post$preds[data_idx,nonchl_cols]) %*% t(post$alpha[,troph_idx, site_idx, nonchl_cols]) # remove chl for autotrophs
    } else {
      lambda <- as.matrix(post$preds[data_idx,]) %*% t(post$alpha[,troph_idx, site_idx, ]) # num per mL
    }
    bmass_conc <- (exp(lambda)*1e6) * rep(mean_bmass/1e9, each = nrow(lambda)) # convert to mgC/m3

    troph_list[[troph]][[site]] <- bmass_conc |> 
      t() |> 
      summarize_pred()

    

    troph_list[[troph]][[site]]$date <- conc_raw$yearmo[data_idx]

  }
  troph_list[[troph]] <- EcotaxaTools::list_to_tib(troph_list[[troph]], 'site')
}

# endregion


all_troph <- EcotaxaTools::list_to_tib(troph_list, 'troph')

# factor role so the shared legend shows all four roles (drop = FALSE), in order
# and with fully written-out names
all_troph$troph <- factor(all_troph$troph, levels = levels(troph_factors))

# regime drawn as contiguous wet/dry run-rectangles below y = 0 (see
# regime_run_bands() in utils.R) rather than a dense daily rug - continuous
# coloured band, no data overlap, tiny file size.
regime_runs <- regime_run_bands(regime)

# shared panel builder. `grp` selects the multi-group community (all roles except
# mixotroph) or the mixotroph column; the mixotroph column drops its y-title and
# subtitle since it shares each row's site with the community panel on its left.
# fill is a single combined scale (roles + wet/dry) so the base regime band and
# the biomass ribbons can coexist; the role legend is supplied separately.
troph_bmass_plot <- function(psite, grp = c('community', 'mixo')) {
  grp <- match.arg(grp)
  sub_data <- if (grp == 'community') {
    all_troph |> filter(site == psite, troph != 'mixotroph')
  } else {
    all_troph |> filter(site == psite, troph == 'mixotroph')
  }

  # regime band is drawn as a strip BELOW y = 0 so it never overlaps the biomass
  # data (which is >= 0); depth is a small fraction of this panel's range, and
  # the y-axis lower limit is extended by the same amount to make room for it
  band_depth <- 0.06 * max(sub_data$high.95, na.rm = TRUE)
  if (!is.finite(band_depth) || band_depth <= 0) band_depth <- 1
  runs <- transform(regime_runs, ymin = -band_depth, ymax = 0)

  ggplot(sub_data) +
    geom_line(aes(x = date, y = median, color = troph)) +
    geom_ribbon(
      aes(x = date, ymin = low.95, ymax = high.95, fill = troph),
      alpha = 0.5
    ) +
    geom_rect(
      data = runs,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = regime),
      inherit.aes = FALSE
    ) +
    labs(
      x = "",
      y = if (grp == 'community') expression('Biomass [mg C m'^-3 * ']') else "",
      subtitle = if (grp == 'community') site_labels[[psite]] else ""
    ) +
    guides(color = 'none', fill = 'none') +
    scale_fill_manual(values = c(troph_cols, regime_cols), drop = FALSE) +
    scale_color_manual(values = troph_cols) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(ylim = c(-band_depth, NA)) +
    theme_pubclean(base_size = 8)
}

# shared role legend, built from all four roles (so every role - including
# Mixotroph - gets a colour swatch) and flattened onto a single row
shared_legend <- get_legend(
  ggplot(all_troph) +
    geom_ribbon(aes(x = date, ymin = low.95, ymax = high.95, fill = troph)) +
    scale_fill_manual(values = troph_cols, labels = troph_labels, drop = FALSE) +
    labs(fill = "") +
    guides(fill = guide_legend(nrow = 1)) +
    theme_pubclean(base_size = 8) +
    theme(legend.position = 'bottom', legend.text = element_text(size = 7))
)

# one complete double-column figure: left column = multi-group community
# (diatoms / dinoflagellates / heterotrophs), right column = mixotrophs, one row
# per site. Flat 2-column ggarrange with the shared role legend at the bottom.
full_plot <- ggarrange(
  troph_bmass_plot('CW', 'community'), troph_bmass_plot('CW', 'mixo'),
  troph_bmass_plot('CE', 'community'), troph_bmass_plot('CE', 'mixo'),
  troph_bmass_plot('AB', 'community'), troph_bmass_plot('AB', 'mixo'),
  troph_bmass_plot('MB', 'community'), troph_bmass_plot('MB', 'mixo'),
  troph_bmass_plot('SC', 'community'), troph_bmass_plot('SC', 'mixo'),
  ncol = 2, nrow = 5,
  align = 'hv',
  common.legend = TRUE,
  legend = 'bottom',
  legend.grob = shared_legend,
  labels = LETTERS
)

ggsave(
  './output/fig07-post_pred_plots.pdf',
  full_plot,
  height = 215, width = 178, units = 'mm',
  dpi = 600
)

