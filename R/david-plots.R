rm(list = ls())
library(ggplot2)


full <- readRDS("./data/02-full_merged.rds")


ggplot(full$conc) +
  geom_point(
    aes(
      x = windspeed,
      y = CHLA_N,
      color = sample_site
    )
  ) +
  geom_smooth(
    aes(
      x = windspeed,
      y = CHLA_N,
      color = sample_site
    ),
    method = 'lm'
  )+
  theme_minimal()

# Diatom biomass vs windspeed 
ggplot(conc_mg |> 
 filter(taxa %in% etx$names$diatom) |> 
 group_by(id, sample_site) |> 
  summarize(
   pgC_L = sum(pgC_L, na.rm = T),
   t = unique(t)
  )
 ) +
  geom_point(
    aes(
      x = t,
      y = log(pgC_L+1),
      color = sample_site
    )
  ) +
  geom_smooth(
    aes(
      x = t,
      y = log(pgC_L+1),
      color = sample_site
    ),
    method = 'lm'
  )+ 
  theme_minimal()
  


####
# NEED TO FIX AND FILTER
###


 
 # Diatom biomass vs salinity 
# diatom_mg = # add filtering sript here

  
 ggplot(diatom_mg) +
   geom_point(
     aes(
       x = sal,
       y = diatom_mg$pgC_L,
       color = sample_site
     )
   ) +
   geom_smooth(
     aes(
       x = sal,
       y = diatom_mg$pgC_L,
       color = sample_site
     )
   )
 
 # Dinoflagellate biomass vs windspeed 
 ggplot(dino_mg) +
   geom_point(
     aes(
       x = windspeed,
       y = dino_mg$pgC_L,
       color = sample_site
     )
   ) +
   geom_smooth(
     aes(
       x = windspeed,
       y = dino_mg$pgC_L,
       color = sample_site
     )
   )
 
  
 # Dinoflagellate biomass vs salinity 
 ggplot(dino_mg) +
   geom_point(
     aes(
       x = sal,
       y = dino_mg$pgC_L,
       color = sample_site
     )
   ) +
   geom_smooth(
     aes(
       x = sal,
       y = dino_mg$pgC_L,
       color = sample_site
     )
   )


#####
# ADD TEMPERUTURE PLOTS
####

###
# Look ciliates
####