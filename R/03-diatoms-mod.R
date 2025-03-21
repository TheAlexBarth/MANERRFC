rm(list= ls())
library(rstan)
library(dplyr)


etx <- readRDS("./data/02-full_merged.RDS")

# trim to diatoms
diat_conc <- etx$conc |> 
  filter(taxa %in% etx$names$diatom)

# flip taxa names 
diat_conc$group <- diat_conc$taxa |> 
  sapply(
    function(x) names(etx$names$diatom)[which(etx$names$diatom == x)]
  ) |> 
  as.factor()

# summarize the counts, numeric concentration, and biomass
# by sample and group (this reduces from 15 taxa to 9 groups)
all_diat <- diat_conc |>
  group_by(sample_id, group, sample_site) |> 
  summarize(
    count = sum(count), 
    num_L = sum(num_L), 
    pgC_L = sum(pgC_L),
    img_vol = unique(img_vol),
    sal = unique(sal),
    temp = unique(t),
    do = unique(DO),
    turb = unique(Turb),
    P = unique(P),
    NH4 = unique(NH4),
    N = unique(N),
    windspeed = unique(windspeed)
  )

par_names <- c('sal', 'temp', 'do','turb','P','NH4','N','windspeed')

# impute to mean
for(par in par_names) {
  if(any(is.na(all_diat[[par]]))) {
    all_diat[[par]][which(is.na(all_diat[[par]]))] <- mean(all_diat[[par]], na.rm = T)
  }
}

diat_indv <- etx$indv |> 
  filter(taxo_name %in% etx$names$diatom)

diat_indv$group <- diat_indv$taxo_name |> 
  sapply(
    function(x) names(etx$names$diatom)[which(etx$names$diatom == x)]
  ) |> 
  as.factor()

# a warning to make sure it is done right
if(any(levels(all_diat$group) != levels(diat_indv$group))) {
  for(i in 1:1e6) {
    print("THE GROUP NAMES ARE WRONG")
  }
}


###################
# MARK: MODEL CONSTRUCTION 
###################

# region \- make predictor ------------

x_counts <- cbind(
  1,
  all_diat$sal,
  all_diat$temp,
  all_diat$do,
  all_diat$turb,
  all_diat$P,
  all_diat$NH4,
  all_diat$N,
  all_diat$windspeed
)

x_indv <- cbind(
  1,
  diat_indv$sal,
  diat_indv$temp,
  diat_indv$P,
  diat_indv$NH4,
  diat_indv$N
)

diat_pois <- list(
  N_obs = nrow(all_diat),
  N_mes = nrow(diat_indv),
  N_groups = length(levels(all_diat$group)),
  n = all_diat$count,
  log_b = log(diat_indv$cmass),
  K_count = ncol(x_counts),
  K_bmass = ncol(x_indv),
  img_vol = all_diat$img_vol,
  
)