rm(list = ls())
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

source('./R/utils.R')
post_role <- readRDS('./data/04-post_micro.RDS')

# region Summarize Effects -----------------------

# amp = sqrt(post_role$alpha[,,,2]^2 + post_role$alpha[,,,3]^2)
alpha_main <- post_role$alpha[,,,-1] # drop int
# alpha_main[,,,1] <- amp # cos to amp


group_levels = list(
  'role' = levels(troph_factors), 
  'site' = levels(site_factors), 
  'w' = names(post_role$preds[,-1])
)

alpha_df <- post_arr_to_df(alpha_main, c('role','site','w'), dim_levels = group_levels)

alpha_sum <- alpha_df |> 
  group_by(role, site, w) |> 
  summarize(
    low = quantile(est, probs = c(0.025)),
    mid = median(est),
    high = quantile(est, probs = c(0.975))
  )


clean_names <- c(
  sin = 'sin_term',
  cos = 'cos_term',
  wind_pca = 'Wind',
  temp = 'Temperature',
  sal = 'Salinity',
  P = 'PO4',
  NH4 = 'NH4',
  N = 'NO23',
  micro = 'Micro',
  nano = 'Nano',
  pico = 'Pico'
)

# endregion -----------------------

# region set NA for some -----------------------

auto_dex <- which(grepl('auto', alpha_sum$role) & alpha_sum$w %in% c('micro','nano','pico'))
alpha_sum$low[auto_dex] <- NA
alpha_sum$mid[auto_dex] <- NA
alpha_sum$high[auto_dex] <- NA

# endregion -----------------------

# region plotting -----------------------

alpha_plot <- function(trrole) {
  sub_data <- alpha_sum |> 
    filter(role == trrole, w != 'sin_term', w != 'cos_term')


  sub_data$w <- clean_names[as.character(sub_data$w)] |> 
    factor(levels = clean_names)

  p = ggplot(sub_data) +
    geom_point(
      aes(x = w, y = mid, color = site),
      position = position_dodge(width = 0.5),
      size = 1.2
    ) +
    geom_segment(
      aes(
        x = w, y = low, yend = high, color = site
      ),
      position = position_dodge(width = 0.5),
      linewidth = 0.3
    ) +
    scale_color_manual(values = site_cols)+
    geom_hline(aes(yintercept = 0), linewidth = 0.3)+
    labs(
      color = "",
      x = "",
      y = "Slope",
      subtitle = trrole
    )+
    theme_pubclean(base_size = 8) + 
    theme(
      axis.title.y = element_text(margin = margin(r = 4)),
      plot.margin = margin(2, 2, 2, 2)
    )
  return(p)
}

diat_auto = alpha_plot('diat_auto')
dino_auto = alpha_plot('dino_auto')
mixotroph = alpha_plot('mixotroph')
grazer = alpha_plot('grazer')



full_plot <- ggarrange(
  diat_auto + theme(axis.text.x = element_text(color = 'transparent')), 
  dino_auto + theme(axis.text.x = element_text(color = 'transparent')),
  mixotroph + theme(axis.text.x = element_text(color = 'transparent')),
  grazer,
  ncol = 1, align = 'v',
  common.legend = TRUE,
  legend = 'bottom'
)
ggsave('./output/fig06-troph_effects.pdf',plot = full_plot, width = 85, height = 120, units = 'mm')


# endregion -----------------------