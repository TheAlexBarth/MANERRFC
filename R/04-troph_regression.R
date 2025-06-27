rm(list = ls())
library(cmdstanr)
library(posterior)
library(dplyr)
library(lubridate)
library(tidyr)
source('./R/utils.R')


etx <- readRDS('./data/02-full_merged.rds')

# region Data prep -----------------------

conc$month <- conc$yearmo |> month()
conc$sin_term <- sin(2*pi*conc$month/12)
conc$cos_term <- sin(2*pi*conc$month/12)

role <- conc$functional_role
site <- conc$sample_site
incld_chl <- ifelse(grepl('_auto',conc$functional_role),0,1)

# endregion -----------------------

conc <- etx$conc

# 