####
# Formatting Ecotaxa Data #####
####


rm(list = ls()) # this wipes your environment
library(EcotaxaTools) # this is my package
library(lubridate)
library(dplyr)
library(ggplot2)
library(tidyr)


## |- Read in data --------------

raw <- read_etx('./data/raw.tsv')
nut <- readRDS('./data/01_swmp_nut.rds')
wq <- readRDS('./data/01_swmp_wq.rds')
wq$sample_site <- wq$StationCode |> substr(4,5) |> toupper()
PDSI <- read.csv('./data/00_PDSI_South-texas.csv')

## |- Prelim formatting ---------

## |-|- Etx format ------------

# filter to all living
names(raw)[which(names(raw) == 'annotation_hierarchy')] <- 'taxo_hierarchy'

living <- raw |> 
  names_keep('living', keep_children = T) 

# renmae to be more simple
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

#change vol ot numeric
living$acq_vol_imaged <- living$acq_vol_imaged |> 
  sapply(function(x) gsub('ml', '',x)) |> 
  as.numeric()

living$acq_dil_fact[is.na(living$acq_dil_fact)] <- 1

## |-|-|- Density Calcs ---------------

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

## |-|- SWMP formatting ---------
# trim to just have dates matching our living periods
wq_summers <- wq |> 
  filter(substr(wq$StationCode, 4,5) %in% c('ce','cw')) |> 
  filter(month %in% c(6:9) & year %in% c(2014:2021))

wq_summer_sum <- wq_summers |> 
  group_by(date = as.Date(wq_summers$DateTimeStamp),
           sample_site) |> 
  summarize(Temp = mean(Temp, na.rm = T),
           Sal = mean(Sal, na.rm = T),
           DO = mean(DO_mgl, na.rm = T),
           Chl = mean(ChlFluor, na.rm = T))


# nut_summers <- nut |> 
#   filter(substr(nut$StationCode, 4,5) %in% c('ce','cw')) |> 
#   filter(month %in% c(6:9) & year %in% c(2014:2021))


###
# Preliminary Plots #############
###

## WQ 02 consumptions

wq$time <- wq$DateTimeStamp |> format('%H:%M:%S')

temp <- wq |>
  filter(month == 6, year == 2021, sample_site == 'CW')

ggplot(temp) +
  geom_line(aes(x = as.POSIXct(time, format = '%H:%M:%S'),
                y = ChlFluor, group = day(DateTimeStamp))) +
  scale_x_datetime(date_breaks = 'hour', date_labels = '%H')

## Data coverage
ggplot(living) +
  geom_point(aes(x = year(date), y = month(date)))+
  facet_grid(~sample_site)


# Salinity over time
wq_summer_sum$moday <- wq_summer_sum$date
year(wq_summer_sum$moday) <- 2000

ggplot(wq_summer_sum) + 
  geom_line(aes(x = moday, y = Sal, color = sample_site)) +
  facet_wrap(~year(wq_summer_sum$date)) + 
  theme_bw() +
  labs(x = '', y = 'Salinity')



# we need to look at the summary data to see visual trends and think how to 

micro_den$moday <- micro_den$date
year(micro_den$moday) <- 2000
micro_trim <- micro_den |> 
  filter(taxa %in% c('Bacillariophyta', 'Dinophyceae', 'Ciliophora'))

ggplot(micro_trim)+
  geom_point(aes(x = moday, y = log(num_L +1), color = taxa, 
                 shape = sample_site)) +
  facet_wrap(~year(micro_trim$date))



## |- Full data set -----------------------
all_data <- micro_den |> 
  left_join(
    wq_summer_sum |> 
      select(Temp, Sal, DO, Chl, date, sample_site),
    by = c('date', 'sample_site')
  )

## \-\- Dist checks ---------
taxa_denPlotter <- function(taxa, remzero = F) {
  data = all_data[all_data$taxa == taxa,]
  
  if(remzero) {
    data = data[data$num_L > 0,]
  }
  
  p <- log(data$num_L + 1) |> 
    density() |> 
    plot()
  
  return(p)
}
  

taxa_denPlotter('Bacillariophyta', F)
taxa_denPlotter('Dinophyceae')
taxa_denPlotter('Ciliophora')

## \-\- Correlation --------



corr_plotter <- function(taxa, col, log = T, remzero = F) {
  
  data = all_data[all_data$taxa == taxa,] |> 
    na.omit()
  
  if(remzero) {
    data = data[data$num_L > 0,]
  }
  
  
  if(log) {
    p <- ggplot(data) +
      geom_point(aes(y = log(num_L + 1), x = data[[col]], color = sample_site)) + 
      labs(y = paste0('Ln(', taxa,')'), x = col)  +
      stat_smooth(aes(y = log(num_L + 1), x = data[[col]]),
                  method = 'lm', se = F) +
      # stat_smooth(aes(y = as.numeric((log(num_L + 1) > 0)), x = data[[col]]),
      #             method = 'glm', method.args = list(family = 'binomial')) +
      theme_bw()
    # print(cor.test(y = log(data[['num_L']] + 1), x = data[[col]], 
    #           method = 'spearman', use = 'complete.obs'))
    # 
    # gam_mod = glm(as.numeric((log(data$num_L +1)>1)) ~ data[[col]],
    #               family = 'binomial')
    # 
    # print(summary(gam_mod))

  } else {
    p <- ggplot(data) +
      geom_point(aes(y = (num_L), x = data[[col]], color = sample_site)) +
      labs(y = taxa, x = col) +
      stat_smooth(aes(y = num_L, x = data[[col]]), 
                  method = 'lm', se = F) +
      theme_bw()
    print(cor.test(y = data[['num_L']], x = data[[col]], 
              method = 'spearman', use = 'complete.obs'))
  }
  return(p)
}

corr_plotter('Ciliophora', 'Chl', log = T, remzero = F)
corr_plotter('Bacillariophyta', 'Sal', log = T, remzero = F)
corr_plotter('Dinophyceae', 'Sal', log = T, remzero = F)

plot(Chl ~ Sal, data = all_data)

###
# groupwise ratio ########
###


## |- Auto vs Hetero -----------------
# 
# hets <- c('Ciliophora', 'Rotifera' ,'heterotroph', 'nauplii', 'Metazoa')
# auts <- c('Bacillariophyta', 'Cyanobacteria', 'Euglenzoa', 'Dinophyceae',
#           "Euglenozoa")
# 
# all_data$aut_het <- NULL
# for(i in 1:nrow(all_data)) {
#   taxon = all_data$taxa[i]
#   if(taxon %in% hets) {
#     all_data$aut_het[i] <- 'het'
#   } else if (taxon %in% auts) {
#     all_data$aut_het[i] <- 'aut'
#   } else if (taxon == 'living') {
#     next
#   } else {
#     stop(paste0('Error at ', i, ' is ', taxon))
#   }
# }

#need ot make a new dataframe
# aut_het_df <- all_data |> 
#   group_by(id, aut_het, sample_site) |> 
#   summarize(num_L = sum(num_L, na.rm = T)) |> 
#   left_join(
#     all_data |> 
#       ungroup() |> 
#       select(id, img_vol, yearmo, date, moday, Temp, Sal, DO, Chl) |> 
#       unique(),
#     by = 'id'
#   )
# 
# # need to pivot
# aut_het_df <- aut_het_df |> 
#   pivot_wider(values_from = num_L, names_from = aut_het)
# 
# 
# aut_het_df$aut_het <- log((aut_het_df$aut + 1) / (aut_het_df$het + 1))
# 
# ggplot(aut_het_df) +
#   geom_point(aes(x = Sal, y = aut_het)) +
#   theme_bw()

## |- Dino vs Diatom -------------------

# need to filter and trim.










