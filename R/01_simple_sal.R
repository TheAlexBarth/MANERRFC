######
# Relationship of Species to Salinity ###
####


rm(list = ls()) # this wipes your environment
library(EcotaxaTools) # this is my package
library(lubridate)
library(dplyr)
library(ggplot2)
library(tidyr)
library(brms)
library(tidybayes)
library(emmeans)

source('./R/utils.R')

## |- Read in data --------------

# will need local adjustment
path = 'C:/Users/David Malcolm/Box/TGCRC Plankton Food Webs/Data/NERRFC'

raw <- read_etx(paste0(path,'/temp_flowcam-full.tsv'))
raw$sample_site <- raw$sample_site |> toupper()
raw <- raw |> filter(sample_id != 'AB_2011-03-08')
raw_CE <- raw |> filter(sample_site == 'CE')

environ <- readRDS('./data/01-environ_clean.rds')
environ$sample_site <- environ$StationCode |> substr(4,5) |> toupper()

# |-|- Etx format ------------

# filter to all living
names(raw)[which(names(raw) == 'annotation_hierarchy')] <- 'taxo_hierarchy'

living <- raw |> 
  names_keep('living', keep_children = T) 

# rename to be more simple
living$taxo_name <- names_to(living, 
                             c('Cyanobacteria', 'Dinophyceae', 'Bacillariophyta',
                               'Euglenozoa','nauplii','Rotifera', 'heterotroph',
                               'Ciliophora', 'Metazoa', 'living'))

## |-|-|- Messy formatting -----------
living$date <- living$sample_date |> 
  sapply(function(x) {
    gsub("-", "", x)
  }) |> 
  as.Date(format = '%Y%m%d')

# get it to be just monthly
living$yearmo <- paste(year(living$date),
                        month(living$date),
                        '01',
                        sep = '-') |> 
  as.Date(format = '%Y-%m-%d')

#change vol to numeric
living$acq_vol_imaged <- living$acq_vol_imaged |> 
  sapply(function(x) gsub('ml', '',x)) |> 
  as.numeric()

living$acq_dil_fact[is.na(living$acq_dil_fact)] <- 1

## |-|-|- Density Calculations ---------------

living$id <- living$sample_id
living$sample_id <- living$acq_id

living_counts <- bin_taxa(living, zooscan = T, force_bins = T)
names(living_counts) <- c('acq_id', 'taxa', 'count')

living_density <- living_counts |> 
  left_join(
    living[, c('acq_id', 'acq_dil_fact', 'id', 'acq_vol_imaged')] |> 
      unique(),
    by = 'acq_id'
  )

living_density$conv_count <- living_density$count / living_density$acq_dil_fact

# final product
micro_den <- living_density |> 
  select(taxa, conv_count, acq_vol_imaged, id) |> 
  group_by(taxa, id) |> 
  summarize(count = sum(conv_count), img_vol = sum(acq_vol_imaged))

micro_den$num_L <- (micro_den$count / micro_den$img_vol) * (1000) #l

micro_den <- micro_den |> 
  left_join(
    living[,c('id', 'yearmo', 'date', 'sample_site')] |> 
      unique(),
    by = 'id'
  )

# \- Join to Environmental

environ_nut_allsite_summers <- environ$nut_avg |> 
  filter(environ$nut_avg$sample_site) %in% c('CE','CW','AB','SC') |> 
  filter(month %in% c(6:9) & year %in% c(2014:2021))

environ_wq_allsite_summers <- environ$wq_sum |> 
  filter(environ$wq_sum$sample_site) %in% c('CE','CW','AB','SC') |> 
  filter(month %in% c(6:9) & year %in% c(2014:2021))

environ_wind_allsite_summers <- environ$wind_avg |> 
  filter(environ$wind_avg$sample_site) %in% c('CE','CW','AB','SC') |> 
  filter(month %in% c(6:9) & year %in% c(2014:2021))

wq_summer_sum <- wq_summers |> 
  group_by(date = as.Date(wq_summers$DateTimeStamp),
           sample_site) |> 
  summarize(Temp = mean(Temp, na.rm = T),
           Sal = mean(Sal, na.rm = T),
           DO = mean(DO_mgl, na.rm = T),
           Chl = mean(ChlFluor, na.rm = T))


## |- Full data set -----------------------
all_data <- micro_den |> 
  left_join(
    wq_summer_sum |> 
      select(Temp, Sal, DO, Chl, date, sample_site),
    by = c('date', 'sample_site')
  )

all_data$sal_scaled <- scale(all_data$Sal)

#####
# Models #########
####


## \- Diatoms --------

just_diat <- all_data |> 
  ungroup() |> 
  filter(taxa == 'Bacillariophyta')

diatom_pois <- brm(
  bf(count ~ sal_scaled + offset(log(img_vol))),
  prior = c(
    set_prior('normal(0,100)', class = 'Intercept'),
    set_prior('normal(0,1)', class = 'b')
  ),
  data = just_diat,
  family = poisson,
  iter = 2000,
  chains = 3,
  thin = 2,
  warmup = 500,
  cores = 10
)

diat_pred <- make_prediction_data(just_diat, 'sal_scaled')

diat_pred$img_vol <- mean(just_diat$img_vol)

diat_pred <- diat_pred |> add_epred_draws(diatom_pois, ndraw = 400, dpar = T)

diat_pred$Sal <- unscale(diat_pred$sal_scaled, mean(just_diat$Sal, na.rm = T), sd(just_diat$Sal, na.rm = T))


ggplot() +
  stat_lineribbon(data = diat_pred,
  aes(x = Sal, y = exp(.epred + log(mean(img_vol)))),
  .width = 0.95, alpha = 0.5) +
  geom_point(data = just_diat,
  aes(x = Sal, y = count)) +
  labs(x = 'Salinity', y = 'Log10(Diatom Count)') +
  scale_y_log10()+
  theme_bw() +
  theme(legend.position = 'none')+
  theme(panel.grid = element_blank(), panel.background = element_blank(),
  axis.title = element_text(size = 8), axis.text = element_text(size = 6))

ggsave('./output/s01_diatom-salinity-pred.pdf', device = 'pdf',
units = 'in', width = 1.5, height = 1.5, dpi= 500)

# \- Ciliates -------------------

just_ciliate <- all_data |> 
  ungroup() |> 
  filter(taxa == 'Ciliophora')

ciliate_pois <- brm(
  bf(count ~ sal_scaled + offset(log(img_vol))),
  prior = c(
    set_prior('normal(0,100)', class = 'Intercept'),
    set_prior('normal(0,1)', class = 'b')
  ),
  data = just_ciliate,
  family = poisson,
  iter = 2000,
  chains = 3,
  thin = 2,
  warmup = 500,
  cores = 10
)

cili_pred <- make_prediction_data(just_ciliate, 'sal_scaled')

cili_pred$img_vol <- mean(just_diat$img_vol)

cili_pred <- cili_pred |> add_epred_draws(ciliate_pois, ndraw = 400, dpar = T)

cili_pred$Sal <- unscale(cili_pred$sal_scaled, mean(just_diat$Sal, na.rm = T), sd(just_diat$Sal, na.rm = T))


ggplot() +
  stat_lineribbon(data = cili_pred,
  aes(x = Sal, y = .epred),
  .width = 0.95, alpha = 0.5) +
  geom_point(data = just_ciliate,
  aes(x = Sal, y = count)) +
  labs(x = 'Salinity', y = 'Log10(Ciliate Count)') +
  scale_y_log10()+
  theme_bw() +
  theme(legend.position = 'none') +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
axis.title = element_text(size = 8), axis.text = element_text(size = 6))
ggsave('./output/s01_ciliate-salinity-pred.pdf', device = 'pdf',
units = 'in', width = 1.5, height = 1.5, dpi= 500)
