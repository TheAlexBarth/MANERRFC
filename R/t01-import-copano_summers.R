rm(list = ls()) # this wipes your environment
library(EcotaxaTools) # this is my package
library(lubridate)
library(dplyr)
library(ggplot2)
library(tidyr)


# region \- Read in data --------------

box_path <- '~/Library/CloudStorage/Box-Box/TGCRC Plankton Food Webs/Data/SWMP Flow Cam Data/Exported/MANERRFC-Copano_summers.tsv'
raw <- read_etx(box_path)

nut <- readRDS('./data/01_swmp_nut.rds')
nut$sample_site <- nut$StationCode |> substr(4,5) |> toupper()
wq <- readRDS('./data/01_swmp_wq.rds')
wq$sample_site <- wq$StationCode |> substr(4,5) |> toupper()

#############################################
# MARK: Etx Formatting ---------------
#############################################



# regioon \- name formatting ---------------------

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

living <- raw |> 
  names_keep(c('living'), keep_children = T)

living$taxo_name <- living |> 
  names_to(c(dinonames, diatom_names, loricated_ciliates, nl_ciliates, 'living'))



# region \-\- Messy formatting -------------------------------------------------------

living$id <- living$sample_id
living$sample_id <- living$acq_id

living$date <- living$sample_date |> 
  sapply(function(x) {
    gsub("-", "", x)
  }) |> 
  as.Date(format = '%Y%m%d')

#change vol ot numeric
living$acq_vol_imaged <- living$acq_vol_imaged |> 
  sapply(function(x) gsub('ml', '',x)) |> 
  as.numeric()

living$acq_dil_fact[is.na(living$acq_dil_fact)] <- 1

# region \-\- Biovolume Calcs ------------------------------------------------------

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



###########################
# MARK: SUMMARIZE ---------------
###########################

living_count <- living |> 
  bin_taxa(zooscan = T, force_bins = T)


names(living_count) <- c('acq_id', 'taxa', 'count')

living_density <- living_count |> 
  left_join(
    living[, c('acq_id', 'acq_dil_fact', 'id', 'acq_vol_imaged')] |> 
      unique(),
    by = 'acq_id'
  )

living_density$conv_count <- living_density$count / living_density$acq_dil_fact

# final counts
micro_den <- living_density |> 
  select(taxa, conv_count, acq_vol_imaged, id) |> 
  group_by(taxa, id) |> 
  summarize(count = sum(conv_count), img_vol = sum(acq_vol_imaged))

micro_den$num_L <- (micro_den$count / micro_den$img_vol) * (1000) #l

micro_den <- micro_den |> 
  left_join(
    living[,c('id', 'date', 'sample_site')] |> 
      unique(),
    by = 'id'
  )



# region \- cmass --------------

living_cmass <- living |> 
  bin_taxa(zooscan = T, func_col = 'cmass', func = sum, force_bins = T)

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

########################################
# MARK: Environmental Format --------
########################################

# averaging values

nut_avg <- nut |> 
  mutate(date = as.Date(DateTimeStamp)) |> 
  select(date, PO4F, NH4F, NO23F, CHLA_N, sample_site) |> 
  group_by(date, sample_site) |> 
  summarize(
    P = mean(PO4F, na.rm = T),
    NH4 = mean(NH4F, na.rm = T),
    N = mean(NO23F, na.rm = T),
    CHLA_N = mean(CHLA_N, na.rm = T)
  ) |> 
  mutate(sample_id = paste(sample_site, date, sep = "_")) |> 
  ungroup()


# water quality daily averages

wq_sum <- wq |> 
  mutate(date = as.Date(DateTimeStamp)) |> 
  group_by(date, sample_site) |> 
  summarize(
          t = mean(Temp, na.rm = T),
          t_max = max(Temp, na.rm = T),
          t_min = min(Temp, na.rm = T),
          s_max = max(Sal, na.rm = T),
          sal = mean(Sal, na.rm = T),
          DO = mean(DO_mgl, na.rm = T),
          DO_min = min(DO_mgl, na.rm = T),
          Chl = mean(ChlFluor, na.rm = T),
          Chl_max = max(ChlFluor, na.rm = T),
          Chl_min = min(ChlFluor, na.rm = T)
          ) |> 
  mutate(sample_id = paste(sample_site, date, sep = '_')) |> 
  ungroup()

# region \- clean up data -------------------

extreme_to_na <- function(vect, sd_away) {
  vect[which(is.infinite(vect))] <- NA
  mv = mean(vect, na.rm = T)
  sv = sd(vect, na.rm = T)

  vout <- vect |> 
    sapply(
      function(x)
      ifelse(x > mv+sd_away*sv, NA, x)
    )
  
  return(vout)
}

# there's some crazy values out there so I'm just removing massive outliers
# that snuck through the QAQC

wq_sum <- wq_sum |> 
  mutate(
    across(c(s_max, sal, DO, DO_min, Chl, Chl_max, Chl_min), ~ extreme_to_na(.,5))
  )



###################################
# MARK: Final Format --------------------
###################################

final_conc <- all_conc |> 
  left_join(
    nut_avg |> 
      select(-c(date, sample_site)),
    by = 'sample_id'
  ) |> 
  left_join(
    wq_sum |> 
      select(-c(date, sample_site)),
    by = 'sample_id'
  )

# this reformats the sample_id to not be acq_id
# doesn't really matter just be aware
final_indv <- living |> 
  as_tibble() |> 
  mutate(sample_id = paste(sample_site, date, sep = "_"))|> 
  select(
    taxo_name,
    taxo_hierarchy,
    pgC = cmass,
    sample_id
  ) |> 
  left_join(
    nut_avg |> 
      select(-c(date, sample_site)),
    by = 'sample_id'
  ) |> 
  left_join(
    wq_sum |> 
      select(-c(date, sample_site)),
    by = 'sample_id'
  )

# \- Save Data -----------------------------------------------

saveRDS(final_conc, './data/t01-copano_concentrations.rds')
saveRDS(final_indv, './data/t01-copano_individals.rds')
saveRDS(
  list(
    dino = dinonames,
    diatom = diatom_names,
    mz = c(loricated_ciliates, nl_ciliates)
  ),
  './data/t01-taxa_names.RDS'
)

# # # curiousity plot:
# # ctemp <- nut_avg |> 
# #   left_join(
# #     wq_sum |> 
# #       select(-date, -sample_site),
# #     by = 'sample_id'
# #   )


# ggplot() +
#   geom_point(
#     data = wq_sum |> 
#       filter(DO < 20 & Chl < 200),
#     aes(
#       x = DO,
#       y = DO_min,
#       color = sample_site
#     ),
#     alpha = 0.5
#   ) +
#   geom_abline(intercept = 0, slope = 1) +
#   geom_hline(yintercept = 3, color = 'red')+
#   geom_vline(xintercept = 3, color = 'red') +
#   theme_minimal()