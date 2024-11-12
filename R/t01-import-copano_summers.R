####
# Formatting Ecotaxa Data #####
####


rm(list = ls()) # this wipes your environment
library(EcotaxaTools) # this is my package
library(lubridate)
library(dplyr)
library(ggplot2)
library(tidyr)


## \- Read in data --------------

box_path <- '~/Library/CloudStorage/Box-Box/TGCRC Plankton Food Webs/Data/SWMP Flow Cam Data/Exported/MANERRFC-Copano_summers.tsv'
raw <- read_etx(box_path)

nut <- readRDS('./data/01_swmp_nut.rds')
wq <- readRDS('./data/01_swmp_wq.rds')



# \- Etx formatting -----------------------------


# \-\- Name Formatting --------------------------------
# filter to all living
names(raw)[which(names(raw) == 'annotation_hierarchy')] <- 'taxo_hierarchy'

# These name vectors should be altered based on final data
# this is good for the summer
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

microzoop_names <- c(
  `Mesodinium` = 'Mesodinium rubrum',
  `Mesodinium` = 'Mesodinium',
  `Other_ciliate` = 'Ciliophora',
  `Other_ciliate` = 't001',
  `Other_ciliate` = 't002',
  `Didinium` = 'Didinium',
  `Spirotrichea_NL` = 'Spirotrichea', #should catch all non-loricated spirotricea
  `Tintinnina` = 'Choreotrichia',
  `Acantharea` = 'Acantharea',
  `Rotifera` = 'Rotifera',
  `Other_metazoan` = 'Metazoa',
  `Other_metazoan` = 'planula',
  `dead_choreo` = 'Choreotrichia X'
)



living <- raw |> 
  names_keep(c('living','t001','t002'), keep_children = T)

living$taxo_name <- living |> 
  names_to(c(dinonames, diatom_names, microzoop_names, 'living', 'othertocheck'))

living$id <- living$sample_id
living$sample_id <- living$acq_id

living$bv <- calc_ellps_vol(living$raw_feret_max, living$raw_feret_min,0.001) #using pixel arg to convert to mm

living_count <- living |> 
  bin_taxa(zooscan = T, force_bins = T)

living_bv <- living |> 
  bin_taxa(zooscan = T, func_col = 'bv', func = sum, force_bins = T)
