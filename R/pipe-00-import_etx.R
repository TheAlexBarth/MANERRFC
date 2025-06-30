####
# Formatting Ecotaxa Data #####
####


############################
# MARK: Initial Prep ---------
############################

rm(list = ls())
library(EcotaxaTools)
library(lubridate)
library(dplyr)
library(ggplot2)
library(tidyr)


# region \- read data ---------------------------
# will need local adjustment
path = '~/Library/CloudStorage/Box-Box/TGCRC Plankton Food Webs/Data/NERRFC'

raw <- read_etx(paste0(path,'/2025_06-revised.tsv'))
raw$sample_site <- raw$sample_site |> toupper()
raw <- raw |> filter(sample_id != 'AB_2011-03-08', abd_diameter >= 20)

taxo_map <- read.csv('./data/00-taxo_map.csv')

######################################################
# MARK: Ecotaxa Format ---------------------------
######################################################

# set to just focus on living
names(raw)[which(names(raw) == 'annotation_hierarchy')] <- 'taxo_hierarchy' # adjust for package

# IMPORTANT NOTE - CAN'T DO THIS
# # IF NOT ALL RAW SAMPLES ARE IN FILTERED LIVING SAMPLE!!!!!
# living2 <- raw |> 
#   names_keep('living', keep_children = T) 
# EASIER TO KEEP not-living and filter later

# rename to be more simple
living <- raw
living$id <- living$sample_id
living$sample_id <- living$acq_id

# region \- name format --------------------------------------
# These name vectors should be altered based on final data
# this is good for the summer but check

living$taxo_name <- living |>
  names_to(c(
    taxo_map$taxo_name,
    'living', 'not-living', 'temporary' #ultimately temporary should be removed
  ))


# final check:
living <- living |> 
  left_join(
    taxo_map |> 
      select(taxo_name, group, functional_role)
  )

living2 <- living |> 
  names_drop(c('living','not-living','temporary'))

table(living2$functional_role, living2$sample_site)


# region \- date format ----------------------------
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

# get it to be just yearly
living$year <- paste(year(living$date), 
                     '01',
                     '01',
                     sep = '-') |>
                       as.Date(format = '%Y-%m-%d')

#date for averaging every month
living$month <- paste('01', 
                      month(living$date),
                      '01',
                      sep = '-') |>
  as.Date(format = '%Y-%m-%d')

#change vol to numeric
living$acq_vol_imaged <- living$acq_vol_imaged |> 
  sapply(function(x) gsub('ml', '',x)) |> 
  as.numeric()

living$acq_dil_fact[is.na(living$acq_dil_fact)] <- 1


# region \- biovolume calcs ----------------------
living$um3 <- (4/3) * pi * (living$abd_diameter/2)^3

living$cmass <- NA
# dinoflagellates
stop("change to base 10 mendenduer")
living$cmass[which(
  living$group == "dinoflagellates"
)] <- exp(-0.353) * living$um3[which(living$group == "dinoflagellates")]^0.864

#little diatoms
living$cmass[which(
  living$group == "diatom" & living$um3 <3000
)] <- exp(-0.541)*living$um3[which(living$group == "diatom" & living$um3 <3000)]^0.811

#big diatoms
living$cmass[which(
  living$group == "diatom" & living$um3 >3000
)] <- exp(-0.933)*living$um3[which(living$group == "diatom" & living$um3 >3000)]^0.881

# loricated ciliates
living$cmass[which(living$group == "lor_ciliate")] <- exp(-0.168)*living$um3[which(living$group == "lor_ciliate")]^0.841

# non-loricated ciliates
living$cmass[which(living$group == "nl_ciliate")] <- exp(-0.639)*living$um3[which(living$group == "nl_ciliate")]^0.984

# NOTE THERE IS A LOT OF NAs LEFT FOR NONFOCUS CATEGORES --WILL BE REMOVED

####################################
# MARK: Ecotaxa Summary --------------
####################################


# region \- count -------------------------
living_count <- living |> 
  bin_taxa(zooscan = T, force_bins = T) |> 
  filter(group %in% taxo_map$taxo_name)

names(living_count) <- c('acq_id', 'taxo_name', 'count')

living_density <- living_count |> 
  left_join(
    living[, c('acq_id', 'acq_dil_fact', 'id', 'acq_vol_imaged')] |> 
      unique(),
    by = 'acq_id'
  )

living_density$conv_count <- living_density$count / living_density$acq_dil_fact

# final product
micro_den <- living_density |> 
  select(taxo_name, conv_count, acq_vol_imaged, id) |> 
  group_by(taxo_name, id) |> 
  summarize(count = sum(conv_count), img_vol = sum(acq_vol_imaged)) |> 
  ungroup() 
micro_den$num_L <- (micro_den$count / micro_den$img_vol) * (1000) #l

micro_den <- micro_den |> 
  left_join(
    living[,c('id', 'yearmo', 'year', 'month', 'date', 'sample_site')] |> 
      unique(),
    by = 'id'
  )|> 
  left_join(
    taxo_map |> 
      select(taxo_name, group, functional_role)
  )
  

# region \- cmass --------------

living_cmass <- living |> 
  bin_taxa(zooscan = T, func_col = 'cmass', func = sum, force_bins = T) |> 
  filter(group %in% taxo_map$taxo_name)


names(living_cmass) <- c('acq_id', 'taxo_name', 'pgC')

cmass_den <- living_cmass |> 
  left_join(
    living[, c('acq_id', 'acq_dil_fact', 'id', 'acq_vol_imaged')] |> 
      unique(),
    by = 'acq_id'
  )

cmass_den$conv_cmass <- cmass_den$pgC / cmass_den$acq_dil_fact

# final counts
cmass_sum <- cmass_den |> 
  select(taxo_name, conv_cmass, acq_vol_imaged, id) |> 
  group_by(taxo_name, id) |> 
  summarize(pgC = sum(conv_cmass), img_vol = sum(acq_vol_imaged)) |> 
  ungroup()

cmass_sum$pgC_L <- (cmass_sum$pgC / cmass_sum$img_vol) * (1000) #l

all_conc <- micro_den |> 
  left_join(
    cmass_sum |> 
      select(id,taxo_name, pgC, pgC_L),
    by = c('id', 'taxo_name')
  ) |> 
  mutate(sample_id = paste(sample_site, date, sep = "_")) |> 
  ungroup() 


saveRDS(
  list(
    conc = all_conc,
    indv = living |> 
      filter(taxo_name %in% taxo_map$taxo_name)
  ),
  "./data/00-ecotaxa_full.rds"
)
