rm(list = ls())
library(ggplot2)
library(tidyr)
library(dplyr)



#############
# MARK: CHOOSE TAXA 
#############

taxa = 'diatom'


##############################
# MARK: LOAD DATA
##############################


full <- readRDS('./data/02-full_merged.rds')
mod = readRDS(paste0('./data/04-mod_output-',taxa,'.RDS'))

source('./R/utils.R')



