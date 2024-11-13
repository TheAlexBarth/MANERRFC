rm(list = ls())
library(EcotaxaTools) # this is my package
library(dplyr)
library(ggplot2)
library(tidyr)
library(brms)
library(tidybayes)
library(emmeans)

conc <- readRDS('./data/t01-copano_concentrations.rds')
name_list = readRDS('./data/t01-taxa_names.rds')

# first just recreating the earlier analyses
## \- Diatoms --------

all_diat <- conc |>
  filter(taxa %in% name_list$diatom) |> 
  group_by(sample_id) |> 
  select(sample_id, count, num_L) |>
  summarize(count = sum(count), num_L = sum(num_L)) |> 
  left_join(
    conc |> 
    select(-c(count, num_L, um3,um3_L, taxa)) |> 
      unique(),
    by = "sample_id"
  )

diatom_pois <- brm(
  bf(count ~ sal + offset(log(img_vol))),
  prior = c(
    set_prior('normal(0,100)', class = 'Intercept'),
    set_prior('normal(0,1)', class = 'b')
  ),
  data = all_diat,
  family = poisson,
  iter = 2000,
  chains = 3,
  thin = 2,
  warmup = 500,
  cores = 10
)

diat_pred <- make_prediction_data(all_diat, 'sal')

diat_pred$img_vol <- mean(all_diat$img_vol)

diat_pred <- diat_pred |> add_epred_draws(diatom_pois, ndraw = 400, dpar = T)

# check as stan:
