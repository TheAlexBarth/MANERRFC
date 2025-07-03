rm(list = ls())
source('./R/utils.R')
post <- readRDS("./data/04-post_micro.RDS")

conc_raw <- readRDS("./data/02-full_merged.rds")$conc
regime <- readRDS('./data/01d-river_regime.rds')
post$alpha |> dim()


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
      lambda <- as.matrix(post$preds[data_idx,1:7]) %*% t(post$alpha[,troph_idx, site_idx, 1:7]) # remove chl for autotrophs
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

all_plotter <- function(psite) {
  p = ggplot(
    all_troph |> 
      filter(site == psite, troph != 'mixotroph')
  ) +
    geom_line(
      aes(
        x = date,
        y = mean,
        color = troph
      )
    )+  
    geom_ribbon(
      aes(
        x = date,
        ymin = low.95,
        ymax = high.95,
        fill = troph
      ),
      alpha = 0.5
    )+
    geom_rug(
      data = regime |> 
        filter(year(yearmo) %in% c(2014:2021)),
      aes(
         x = Date,
         color = regime
      )
    )+
    labs(x = "", y = 'Biomass [mgC/m3]', fill = "", subtitle = psite)+
    guides(color = 'none')+
    scale_fill_manual(values = troph_cols)+
    scale_color_manual(values = c(troph_cols, regime_cols)) +
    theme_pubclean(base_size = 8)
  return(p)
}


cw_plot <- all_plotter('CW')
ce_plot <- all_plotter("CE")
ab_plot <- all_plotter("AB")
mb_plot <- all_plotter('MB')
sc_plot <- all_plotter("SC")

full_plot <- ggarrange(
  cw_plot,
  ce_plot,
  ab_plot,
  mb_plot,
  sc_plot,
  ncol = 1,
  align = 'v',
  common.legend = T
)

ggsave(
  './output/fig07-post_pred_plots.pdf',
  full_plot,
  height = 270, width = 170, units = 'mm',
  dpi= 600
)


# region Mixoplot -----------------------
mixo_plotter <- function(psite) {
  p = ggplot(
    all_troph |> 
      filter(site == psite, troph == 'mixotroph')
  ) +
    geom_line(
      aes(
        x = date,
        y = mean,
        color = troph
      )
    )+  
    geom_ribbon(
      aes(
        x = date,
        ymin = low.95,
        ymax = high.95,
        fill = troph
      ),
      alpha = 0.5
    )+
    geom_rug(
      data = regime |> 
        filter(year(yearmo) %in% c(2014:2021)),
      aes(
         x = Date,
         color = regime
      )
    )+
    labs(x = "", y = 'Biomass [mgC/m3]', fill = "", subtitle = psite)+
    guides(color = 'none')+
    scale_fill_manual(values = troph_cols)+
    scale_color_manual(values = c(troph_cols, regime_cols)) +
    theme_pubclean(base_size = 8)
  return(p)
}


cw_mplot <- mixo_plotter('CW')
ce_mplot <- mixo_plotter("CE")
ab_mplot <- mixo_plotter("AB")
mb_mplot <- mixo_plotter('MB')
sc_mplot <- mixo_plotter("SC")

full_mplot <- ggarrange(
  cw_mplot,
  ce_mplot,
  ab_mplot,
  mb_mplot,
  sc_mplot,
  ncol = 1,
  align = 'v',
  common.legend = T
)

ggsave(
  './output/fig07b-post_pred_plots.pdf',
  full_mplot,
  height = 270, width = 170, units = 'mm',
  dpi= 600
)


# endregion -----------------------

