#####
# Example comparison wet dry
#######

rm(list = ls())
library(lubridate)
library(dplyr)
# do diatoms change in wet vs dry periods?

etx <- readRDS('./data/02-full_merged.rds')
env <- readRDS('./data/01-environ_clean.RDS')


env$PDSI$yearmo <- as.Date(env$PDSI$Date)

## QUICK CREATE STATE:

# choose to set based on pdsi index
# env$PDSI$state = 'normal'
# env$PDSI$state[env$PDSI$pdsi <= -2] <- 'dry' # index based on psdi and assign state.
# env$PDSI$state[env$PDSI$pdsi > 0] <- 'wet'
# add turbidity as well


plank_pdsi <- etx$conc |> 
  filter(
    taxa %in% etx$names$diatom
  ) |> 
  group_by(
    sample_site, yearmo
  ) |> 
  summarize(
    mean_c = sum(pgC_L),
    sd_c = sd(pgC_L)
  ) |> 
  left_join(
    env$PDSI |> 
      select(
        yearmo,
        pdsi
      )
  )


ggplot() +
  geom_point(
    data = plank_pdsi,
    aes(
      x = pdsi,
      y = mean_c,
      color = sample_site
    )
  )
