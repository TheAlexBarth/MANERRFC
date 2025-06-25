################
# Figures for relative proportion
################


rm(list = ls())

library(ggplot2)
library(dplyr)


# region \- data prep -----------------

etx <- readRDS("./data/02-full_merged.RDS")


# region \- Rel by site ----------------

total_by_site <- etx$conc |> 
  group_by(functional_role, yearmo, sample_site) |> 
  summarize(
    pgC_L = sum(pgC_L),
    num_L = sum(num_L)
  ) |> 
  left_join(
    dm_total
  )

site_group_plot <- function(site, metric) {
  sub_data <- total_by_site |> 
    filter(sample_site == site)
  p = ggplot() +
    geom_bar(
      aes(
        x = sub_data$yearmo, y = sub_data[[metric]],
        fill = sub_data$functional_role
      ),
      position = 'stack', stat = 'identity'
    ) +
    labs(
      x = "", y= "", fill = "", 
      subtitle = paste0(subtitle = paste0(site, " - ", metric))
    )+
    scale_y_continuous(limits = c(0,2.1e6))+
    theme_minimal()
}

# # biomass
# site_biomass_list <- list()
# for(site in total_by_site$sample_site) {
#   site_biomass_list[[site]] <- site_group_plot(site, 'pgC_L')
# }
# site_biomass_list$CW
# site_biomass_list$CE
# site_biomass_list$AB
# site_biomass_list$MB
# site_biomass_list$SC

# numeric conc
site_conc_list <- list()
for(site in total_by_site$sample_site) {
  site_conc_list[[site]] <- site_group_plot(site, 'num_L')
}
site_conc_list$CW
site_conc_list$CE
site_conc_list$AB
site_conc_list$MB
site_conc_list$SC

# endregion


# region \- seasonal ----------
season_by_site <- total_by_site |> 
  group_by(
    sample_site, functional_role, 
    month = month(yearmo), drought
  ) |> 
  summarize(
    mean_den = mean(num_L),
    sd_den = sd(num_L)
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
      position = 'stack', stat = 'identity'
    ) +
    labs(
      x = "", y= "", fill = "", 
      subtitle = paste0(subtitle = site)
    )+
    facet_wrap(~ drought) +
    theme_minimal()
}

cw_site <- season_site('CW')
ce_site <- season_site("CE")
ab_site <- season_site('AB')
mb_site <- season_site('MB')
sc_site <- season_site("SC")

cw_site
ce_site
ab_site
mb_site
sc_site




# endregion

# region \- Supplemental Nutrient Plots -------------------


# endregion