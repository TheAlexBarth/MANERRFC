################
# Figures for relative proportion
################


rm(list = ls())

source('./R/utils.R')
library(ggplot2)
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
    labs(
      x = "", y= "", fill = "", 
      subtitle = site
    )+
    geom_rug(
      data = regime |> 
        filter(year(Date) %in% c(min(year(sub_data$yearmo)):max(year(sub_data$yearmo)))),
      aes(
        x = Date,
        color = regime
      )
    )+
    scale_fill_manual(values = troph_cols) +
    scale_color_manual(values = regime_cols) +
    scale_y_continuous(limits = c(0, 1500))+
    labs(x = "", y = "Number per mL", fill = "")+
    guides(color = 'none')+
    theme_minimal()
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
full_times, width = 170, height = 250, units = "mm",
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