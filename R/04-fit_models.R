########
#  This script is run after the 03 script which was used to select the best model
#  This will loop and save the individual models
########


rm(list = ls())
source('./R/utils-mod_making.R')

etx <- readRDS("./data/02-full_merged.RDS")

# region \- mz ------------------------
mz_count_preds <- c('sal','t')
fit_model(etx, 'mz',mz_count_preds)

# region \- dinos -------------
dino_count_preds <- c('sal','t','do','turb', 'windspeed')
fit_model(etx, 'dino',dino_count_preds)

# region \- diatoms ---------
diatom_count_preds <-  c('sal','t','turb', 'windspeed')
fit_model(etx, 'diatom', diatom_count_preds)

