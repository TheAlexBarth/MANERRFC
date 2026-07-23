rm(list = ls())
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

source('./R/utils.R')
post_chl <- readRDS('./data/03-post_chl.RDS')


# region  Summarize effects -----------------------


# amp = sqrt(post_chl$beta[,,,2]^2 + post_chl$beta[,,,3]^2)
beta_main <- post_chl$beta[,,,-c(1,2,3)] # drop sin term and intecept and cos
# beta_main[,,,1] <- amp # replace cos with amp



group_levels = list(
  'frac' = levels(size_factors), 
  'site' = levels(site_factors), 
  'x' = c(names(post_chl$data$X_scaled[,-c(1,2,3)]))
)

beta_df <- post_arr_to_df(beta_main, c('frac','site','x'), dim_levels = group_levels)


beta_sum <- beta_df |> 
  group_by(frac, site, x) |> 
  summarize(
    low = quantile(est, probs = c(0.025)),
    mid = median(est),
    high = quantile(est, probs = c(0.975))
  )


clean_names <- c(
  seasonal = 'Season',
  wind = 'Wind',
  temp = 'Temperature',
  sal = 'Salinity',
  P = 'PO4',
  NH4 = 'NH4',
  N = 'NO23',
  SiOH = 'SiO4'
)

beta_plot <- function(size_frac) {
  sub_data <- beta_sum |> 
    filter(frac == size_frac)

  sub_data$x <- clean_names[as.character(sub_data$x)] |> 
    factor(levels = clean_names)

  p = ggplot(sub_data) +
    geom_point(
      aes(x = x, y = mid, color = site),
      position = position_dodge(width = 0.5),
      size = 1.2
    ) +
    geom_segment(
      aes(
        x = x, y = low, yend = high, color = site
      ),
      position = position_dodge(width = 0.5),
      linewidth = 0.3
    ) +
    scale_color_manual(values = site_cols)+  # abbreviated station names in legend
    geom_hline(aes(yintercept = 0), linewidth = 0.3)+
    labs(
      color = "",
      x = "",
      y = "Slope",
      subtitle = size_labels[[size_frac]]
    )+
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 1.5)))+
    theme_pubclean(base_size = 8) +
    theme(
      axis.title.y = element_text(margin = margin(r = 4)),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
      legend.text = element_text(size = 6.5),
      legend.key.size = unit(9, 'pt'),
      plot.margin = margin(2, 2, 2, 2)
    ) +
    scale_y_continuous(limits = c(-2.12,2))
  return(p)
}

pico = beta_plot('pico')
nano = beta_plot('nano')
micro = beta_plot('micro')


# endregion -----------------------


# only the bottom (pico) panel shows the x-axis; fully blank it on the upper two
# (not just transparent text) so the angled labels don't reserve empty space and
# squeeze the panels
hide_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

full_plot <- ggarrange(
  micro + hide_x,
  nano + hide_x,
  pico,
  ncol = 1, align = 'v',
  common.legend = TRUE,
  legend = 'bottom',
  labels = LETTERS
)
ggsave('./output/fig04-chl_effects.pdf', plot = full_plot,
  width = 85, height = 130, units = 'mm', dpi = 600)
