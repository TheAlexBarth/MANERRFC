rm(list = ls())
library(EcotaxaTools) # this is my package
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(rstan)
library(lubridate)
set.seed(0213)
source('./R/utils.R')


# read conc data
conc <- readRDS('./data/t01-copano_concentrations.rds')
name_list = readRDS('./data/t01-taxa_names.rds')

# read individual data
indv <- readRDS('./data/t01-copano_individals.rds')


#
########################
# MARK: Diatoms --------
########################


# trim to diatoms
diat_conc <- conc |> 
  filter(taxa %in% name_list$diatom)

# flip taxa names 
diat_conc$group <- diat_conc$taxa |> 
  sapply(
    function(x) names(name_list$diatom)[which(name_list$diatom == x)]
  ) |> 
  as.factor()


all_diat <- diat_conc |>
  group_by(sample_id, group) |> 
  summarize(count = sum(count), num_L = sum(num_L), pgC_L = sum(pgC_L)) |> 
  left_join(
    diat_conc |> 
    select(-c(count, num_L, pgC,pgC_L, taxa, group)) |> 
      unique(),
    by = "sample_id"
  )

# impute to mean
all_diat$sal[which(is.na(all_diat$sal))] <- mean(all_diat$sal, na.rm = T)

diat_indv <- indv |> 
  filter(taxo_name %in% name_list$diatom)

diat_indv$group <- diat_indv$taxo_name |> 
  sapply(
    function(x) names(name_list$diatom)[which(name_list$diatom == x)]
  ) |> 
  as.factor()


# region \- stan fit ----------
diat_indv$sal[is.na(diat_indv$sal)] <- mean(diat_indv$sal, na.rm = T)

# don't use mix-multiple since there is no effect of temperature
diat_pois <- list(
  N_obs = nrow(all_diat),
  N_mes = nrow(diat_indv),
  N_groups = length(unique(all_diat$group)),
  n = all_diat$count,
  group = as.numeric(all_diat$group),
  img_vol = all_diat$img_vol/1000, # convert to L,
  sal = all_diat$sal,
  group_counts = as.numeric(all_diat$group),
  group_bio = as.numeric(diat_indv$group),
  log_b = log(diat_indv$pgC),
  sal_b = diat_indv$sal,
  temp_b = diat_indv$t,
  obs_biomass_conc = all_diat$pgC_L
)

diat_fit <- stan(
  file = './stan/03b-mix-sal_only-hierachical.stan',
  data = diat_pois,
  chains = 4, iter = 5000, warmup = 500, cores = 14
)
# region \- save -------------
saveRDS(
  list(mod = diat_fit, data = all_diat, indv = diat_indv),
   "./data/msi-diat_fit.rds"
)


##########################################
# MARK: MZ ---------------------------
##########################################



# trim to mzoms
mz_conc <- conc |> 
  filter(taxa %in% name_list$mz)

# flip taxa names 
mz_conc$group <- mz_conc$taxa |> 
  sapply(
    function(x) names(name_list$mz)[which(name_list$mz == x)]
  ) |> 
  as.factor()


all_mz <- mz_conc |>
  group_by(sample_id, group) |> 
  summarize(count = sum(count), num_L = sum(num_L), pgC_L = sum(pgC_L)) |> 
  left_join(
    mz_conc |> 
    select(-c(count, num_L, pgC, pgC_L, taxa, group)) |> 
      unique(),
    by = "sample_id"
  )


mz_indv <- indv |> 
  filter(taxo_name %in% name_list$mz)

mz_indv$group <- mz_indv$taxo_name |> 
  sapply(
    function(x) names(name_list$mz)[which(name_list$mz == x)]
  ) |> 
  as.factor()
  


# impute to mean
all_mz$sal[which(is.na(all_mz$sal))] <- mean(all_mz$sal, na.rm = T)
mz_indv$t[which(is.na(mz_indv$t))] <- mean(mz_indv$t, na.rm = T)


# region \- stan fit ----------
mz_indv$sal[is.na(mz_indv$sal)] <- mean(mz_indv$sal, na.rm = T)

# don't use mix-multiple since there is no effect of temperature
mz_pois <- list(
  N_obs = nrow(all_mz),
  N_mes = nrow(mz_indv),
  N_groups = length(unique(all_mz$group)),
  n = all_mz$count,
  group = as.numeric(all_mz$group),
  img_vol = all_mz$img_vol/1000, # convert to L,
  sal = all_mz$sal,
  group_counts = as.numeric(all_mz$group),
  group_bio = as.numeric(mz_indv$group),
  log_b = log(mz_indv$pgC),
  sal_b = mz_indv$sal,
  temp_b = mz_indv$t,
  obs_biomass_conc = all_mz$pgC_L
)

mz_fit <- stan(
  file = './stan/03b-mix-sal_only-hierachical.stan',
  data = mz_pois,
  chains = 4, iter = 5000, warmup = 500, cores = 14
)
# region \- save -------------
saveRDS(
  list(
    mod = mz_fit,
    data = all_mz,
    indv = mz_indv
  ),
  './data/msi-mz_fit.rds'
)

################################
# MARK: Dinos ---------------------------
##########################################

# trim to dinooms
dino_conc <- conc |> 
  filter(taxa %in% name_list$dino)

# flip taxa names 
dino_conc$group <- dino_conc$taxa |> 
  sapply(
    function(x) names(name_list$dino)[which(name_list$dino == x)]
  ) |> 
  as.factor()


all_dino <- dino_conc |>
  group_by(sample_id, group) |> 
  summarize(count = sum(count), pgC = sum(num_L), pgC_L = sum(pgC_L)) |> 
  left_join(
    dino_conc |> 
    select(-c(count, num_L, pgC,pgC_L, taxa, group)) |> 
      unique(),
    by = "sample_id"
  )


dino_indv <- indv |> 
  filter(taxo_name %in% name_list$dino)

dino_indv$group <- dino_indv$taxo_name |> 
  sapply(
    function(x) names(name_list$dino)[which(name_list$dino == x)]
  ) |> 
  as.factor()
  


# impute to mean
all_dino$sal[which(is.na(all_dino$sal))] <- mean(all_dino$sal, na.rm = T)
dino_indv$t[which(is.na(dino_indv$t))] <- mean(dino_indv$t, na.rm = T)


# region \- stan fit ----------
dino_indv$sal[is.na(dino_indv$sal)] <- mean(dino_indv$sal, na.rm = T)

# don't use mix-multiple since there is no effect of temperature
dino_pois <- list(
  N_obs = nrow(all_dino),
  N_mes = nrow(dino_indv),
  N_groups = length(unique(all_dino$group)),
  n = all_dino$count,
  group = as.numeric(all_dino$group),
  img_vol = all_dino$img_vol/1000, # convert to L,
  sal = all_dino$sal,
  group_counts = as.numeric(all_dino$group),
  group_bio = as.numeric(dino_indv$group),
  log_b = log(dino_indv$pgC),
  sal_b = dino_indv$sal,
  temp_b = dino_indv$t,
  obs_biomass_conc = all_dino$pgC_L
)

dino_fit <- stan(
  file = './stan/03b-mix-sal_only-hierachical.stan',
  data = dino_pois,
  chains = 4, iter = 5000, warmup = 500, cores = 14
)


# region \- save -------------
saveRDS(
  list(
    mod = dino_fit,
    data = all_dino,
    indv = dino_indv
  ),
  './data/msi-dino_fit.rds'
)