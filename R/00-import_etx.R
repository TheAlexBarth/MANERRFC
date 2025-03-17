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
 path = 'C:/Users/David Malcolm/Box/TGCRC Plankton Food Webs/Data/NERRFC'
#path = '~/Library/CloudStorage/Box-Box/TGCRC Plankton Food Webs/Data/NERRFC'

raw <- read_etx(paste0(path,'/temp_flowcam-full.tsv'))
raw$sample_site <- raw$sample_site |> toupper()
raw <- raw |> filter(sample_id != 'AB_2011-03-08')

######################################################
# MARK: Ecotaxa Format ---------------------------
######################################################

# set to just focus on living
names(raw)[which(names(raw) == 'annotation_hierarchy')] <- 'taxo_hierarchy' # adjust for package

# IMPORTANT NOTE - CAN'T DO THIS
# IF NOT ALL RAW SAMPLES ARE IN FILTERED LIVING SAMPLE!!!!!
# living <- raw |> 
#   names_keep('living', keep_children = T) 
# EASIER TO KEEP not-living and filter later

# rename to be more simple
living <- raw
living$id <- living$sample_id
living$sample_id <- living$acq_id

# region \- name format --------------------------------------
# These name vectors should be altered based on final data
# this is good for the summer but check
# table(raw$taxo_name) |> sort()

dinonames <- c(
  `Other_dino` = 'Dinophyceae',
  `Dinophysis` = 'Dinophysis',
  `Other_armored` = 'Gonyaulacales',
  `Other_unarmored` = 'Gymnodinales',
  `Akashiwo` = 'Akashiwo',
  `Margalef` = 'Cochlodinium',
  `Gymnodinium` = 'Gymnodinium',
  `Gyrodinium` = 'Gyrodinium',
  `Karenia` = 'Karenia',
  `Other_armored` = 'Peridinales',
  `Other_unarmored` = 'Phalacroma',
  `Prorocentrum` = 'Prorocentrum'
)

diatom_names <- c(
  `solo_pennate` = 'Bacillariophyta',
  `Asterionellopsis` = 'Asterionellopsis',
  `Baxillaria` = 'Bacillaria paxillifera',
  `solo_pennate` = 'Cylindrotheca',
  `solo_pennate` = 'Entomoneis',
  `solo_pennate` = 'Pleurosigma',
  `Thalassionema` = 'Thalassionema',
  `solo_pennate` = 'tropidoneis',
  `centric_chain` = 'Lioloma',
  `centric_chain` = 'Thalassiosira',
  `centric_chain` = 'centric_chain',
  `solo_centric` = 'centric',
  `pennate_chain` = 'Pseudo-nitzschia',
  `pennate_chain` = 'Striatella',
  `Mediophyceae` = 'Mediophyceae',
  `Coscinodiscids` = 'Coscinodiscophytina'
)

loricated_ciliates <- c(
  `Tintinnina` = 'Choreotrichia',
  `Tintinnina` = 'Choreotrichia X'
)

nl_ciliates <- c(
  `Mesodinium` = 'Mesodinium rubrum',
  `Mesodinium` = 'Mesodinium',
  `Other_ciliate` = 'Ciliophora',
  `Other_ciliate` = 't001',
  `Other_ciliate` = 't002',
  `Didinium` = 'Didinium',
  `Spirotrichea_NL` = 'Spirotrichea' #should catch all non-loricated spirotricea
)

living$taxo_name <- living |>
  names_to(c(
    dinonames, diatom_names, loricated_ciliates, nl_ciliates,
    'living', 'not-living', 'temporary' #ultimately temporary should be removed
  ))

# final check:
 table(living$taxo_name) |> sort()

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

#change vol to numeric
living$acq_vol_imaged <- living$acq_vol_imaged |> 
  sapply(function(x) gsub('ml', '',x)) |> 
  as.numeric()

living$acq_dil_fact[is.na(living$acq_dil_fact)] <- 1


# region \- biovolume calcs ----------------------
living$um3 <- (4/3) * pi * (living$abd_diameter/2)^3

living$cmass <- NA
# dinoflagellates
living$cmass[which(
  living$taxo_name %in% dinonames
)] <- exp(-0.353) * living$um3[which(living$taxo_name %in% dinonames)]^0.864

#little diatoms
living$cmass[which(
  living$taxo_name %in% diatom_names & living$um3 <3000
)] <- exp(-0.541)*living$um3[which(living$taxo_name %in% diatom_names & living$um3 <3000)]^0.811

#big diatoms
living$cmass[which(
  living$taxo_name %in% diatom_names & living$um3 >3000
)] <- exp(-0.933)*living$um3[which(living$taxo_name %in% diatom_names & living$um3 >3000)]^0.881

# loricated ciliates
living$cmass[which(living$taxo_name %in% loricated_ciliates)] <- exp(-0.168)*living$um3[which(living$taxo_name %in% loricated_ciliates)]^0.841

# non-loricated ciliates
living$cmass[which(living$taxo_name %in% nl_ciliates)] <- exp(-0.639)*living$um3[which(living$taxo_name %in% nl_ciliates)]^0.984

# NOTE THERE IS A LOT OF NAs LEFT FOR NONFOCUS CATEGORES --WILL BE REMOVED

####################################
# MARK: Ecotaxa Summary --------------
####################################


# region \- count -------------------------
living_count <- living |> 
  bin_taxa(zooscan = T, force_bins = T) |> 
  filter(group %in% c(dinonames, diatom_names, loricated_ciliates, nl_ciliates))

names(living_count) <- c('acq_id', 'taxa', 'count')

living_density <- living_count |> 
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

# region \- cmass --------------

living_cmass <- living |> 
  bin_taxa(zooscan = T, func_col = 'cmass', func = sum, force_bins = T) |> 
  filter(group %in% c(dinonames, diatom_names, loricated_ciliates, nl_ciliates))


names(living_cmass) <- c('acq_id', 'taxa', 'pgC')

cmass_den <- living_cmass |> 
  left_join(
    living[, c('acq_id', 'acq_dil_fact', 'id', 'acq_vol_imaged')] |> 
      unique(),
    by = 'acq_id'
  )

cmass_den$conv_cmass <- cmass_den$pgC / cmass_den$acq_dil_fact

# final counts
cmass_sum <- cmass_den |> 
  select(taxa, conv_cmass, acq_vol_imaged, id) |> 
  group_by(taxa, id) |> 
  summarize(pgC = sum(conv_cmass), img_vol = sum(acq_vol_imaged))

cmass_sum$pgC_L <- (cmass_sum$pgC / cmass_sum$img_vol) * (1000) #l

all_conc <- micro_den |> 
  left_join(
    cmass_sum |> 
      select(id,taxa, pgC, pgC_L),
    by = c('id', 'taxa')
  ) |> 
  mutate(sample_id = paste(sample_site, date, sep = "_")) |> 
  ungroup()


saveRDS(
  list(
    names = list(
      dino = dinonames,
      diatom = diatom_names,
      mz = c(loricated_ciliates, nl_ciliates)
    ),
    conc = all_conc,
    indv = living |> 
      filter(taxo_name %in% c(dinonames, diatom_names, loricated_ciliates, nl_ciliates))
  ),
  "./data/00-ecotaxa_full.rds"
)
