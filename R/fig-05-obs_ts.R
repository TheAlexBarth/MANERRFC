################
# Figures for relative proportion
################


rm(list = ls())

source('./R/utils.R')
library(ggplot2)
library(ggpubr)
library(dplyr)
library(lubridate)


# region \- data prep -----------------

etx <- readRDS("./data/02-full_merged.rds")
regime <- readRDS("./data/01d-river_regime.RDS")


# region \- Rel by site ----------------

total_by_site <- etx$conc |>
  group_by(functional_role, id, sample_site, yearmo) |>
  summarize(
    num_L = sum(count) / unique(img_vol)
  ) |>
  ungroup() |>
  group_by(functional_role, sample_site, yearmo) |>
  summarize(num_mL = mean(num_L))

# wet/dry regime drawn as a continuous band below y = 0 (contiguous run
# rectangles, see regime_run_bands() in utils.R) instead of a dense daily rug
ymax_num <- 1500
regime_depth <- 0.06 * ymax_num
regime_runs <- transform(regime_run_bands(regime), ymin = -regime_depth, ymax = 0)

site_group_plot <- function(site) {
  sub_data <- total_by_site |>
    filter(sample_site == site)
  p = ggplot(sub_data) +
    geom_bar(
      aes(
        x = yearmo, y = num_mL,
        fill = functional_role
      ),
      position = 'stack', stat = 'identity'
    ) +
    geom_rect(data = subset(regime_runs, regime == 'Dry'),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = regime_cols[['Dry']], inherit.aes = FALSE) +
    geom_rect(data = subset(regime_runs, regime == 'Wet'),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = regime_cols[['Wet']], inherit.aes = FALSE) +
    scale_fill_manual(values = troph_cols, labels = troph_labels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = scales::label_comma())+
    coord_cartesian(ylim = c(-regime_depth, ymax_num))+
    guides(fill = guide_legend(nrow = 2))+
    labs(
      x = "", y = expression("Number ["*mL^-1*"]"), fill = "",
      subtitle = site_labels[[site]]
    )+
    theme_minimal(base_size = 8) +
    theme(
      legend.text = element_text(size = 6.5),
      legend.key.size = unit(9, 'pt'),
      plot.subtitle = element_text(size = 8)
    )
}


# numeric conc
site_conc_list <- list()
for(site in total_by_site$sample_site) {
  site_conc_list[[site]] <- site_group_plot(site)
}

full_times <- ggarrange(
  site_conc_list[[1]],
  site_conc_list[[2]],
  site_conc_list[[3]],
  site_conc_list[[4]],
  site_conc_list[[5]],
  common.legend = TRUE,
  ncol = 1,
  align = 'v'
)

ggsave('./output/fig05-obs_num_conc.pdf',
full_times, width = 89, height = 220, units = "mm",
dpi = 600
)



# endregion


# region \- seasonal ----------


season_by_site <- total_by_site |> 
  left_join(
    regime |> select(yearmo, regime),
    relationship = 'many-to-many'
  ) |> 
  group_by(
    sample_site, functional_role, 
    month = month(yearmo), regime
  ) |> 
  summarize(
    mean_den = mean(num_mL),
    sd_den = sd(num_mL)
  )


season_site <- function(site) {
  sub_data <- season_by_site |> 
    filter(sample_site == site)
  p = ggplot(sub_data) +
    geom_bar(
      aes(
        x = month, y = mean_den,
        fill = functional_role
      ),
      position = position_dodge(0.5), stat = 'identity'
    ) +
    geom_errorbar(
      aes(
        x = month, ymax = mean_den + sd_den,
        ymin = mean_den,
        color = functional_role
      ),
      position = position_dodge(0.5)
    )+
    labs(
      x = "", y= "", fill = "", 
      subtitle = paste0(subtitle = site)
    )+
    facet_wrap(~ regime) +
    scale_color_manual(values = troph_cols) +
    scale_fill_manual(values = troph_cols) +
    guides(color = "none")+
    labs(x = 'Month')+
    theme_minimal()
}

cw_site <- season_site('CW')
ce_site <- season_site("CE")
ab_site <- season_site('AB')
mb_site <- season_site('MB')
sc_site <- season_site("SC")

seasonal_time <- ggarrange(
  cw_site, ce_site, ab_site, mb_site, sc_site,
  ncol = 1,
  common.legend = T, align = 'v',
  labels = LETTERS
)
ggsave(
  './output/s-seasonal_numeric_mean.pdf',seasonal_time,
  width = 170, height = 280, units = 'mm', dpi = 600
)


# endregion