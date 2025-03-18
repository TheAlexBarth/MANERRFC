#####
# Example comparison wet dry
#######

rm(list = ls())
library(lubridate)
library(dplyr)
library(ggplot2)

# do diatoms change in wet vs dry periods?

etx <- readRDS('./data/02-full_merged.rds')
env <- readRDS('./data/01-environ_clean.RDS')


env$PDSI$yearmo <- as.Date(env$PDSI$Date)

## QUICK CREATE STATE:

# choose to set based on pdsi index
env$PDSI$state = 'normal'
env$PDSI$state[env$PDSI$pdsi <= -2] <- 'dry' # index based on pdsi and assign state.
env$PDSI$state[env$PDSI$pdsi > 0] <- 'wet'
# add turbidity as well

### All Diatoms Biomass Vs PDSI
diatom_pdsi_all <- etx$conc |> 
  filter(
    taxa %in% etx$names$diatom
  ) |> 
  group_by(
    sample_site, yearmo
  ) |> 
  summarize(
    mean_c = sum(pgC_L, na.rm = T),
    sd_c = sd(pgC_L, na.rm = T)
  ) |> 
  left_join(
    env$PDSI |> 
      select(
        yearmo,
        pdsi
      )
  )

### Dry Diatoms Biomass vs PDSI
diatom_pdsi_dry <- etx$conc |> 
  filter(
    taxa %in% etx$names$diatom
  ) |> 
  group_by(
    sample_site, yearmo
  ) |> 
  summarize(
    mean_c = sum(pgC_L, na.rm = T),
    sd_c = sd(pgC_L, na.rm = T)
  ) |> 
  left_join(
    env$PDSI |> 
      filter(state == 'dry') |> #Filters only 'dry' PDSI values
      select(
        yearmo,
        pdsi,
        windspeed
      )
  )


### Plot All Diatom Biomass vs PDSI
ggplot() +
  geom_point(
    data = diatom_pdsi_all,
    aes(
      x = pdsi,
      y = log(mean_c),
      color = sample_site
    )
  ) +
  geom_smooth(
    data = diatom_pdsi_all,
    aes(
      x = pdsi,
      y = log(mean_c),
      color = sample_site
    ),
    method = 'lm'
  )+
  theme_minimal()
