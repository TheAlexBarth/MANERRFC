######
# Environmental Formatting
#
# Imports SWMP water-quality trip records and formats monthly site means for
# temp, salinity, PO4 (P), NH4, NO23 (N), CHLA, and SiOH (silicate). SiOH was
# added per reviewer request; it is folded in here (formerly pipe-01e) because
# it comes from the same raw files and QAQC conventions. Unlike temp/sal it has
# near-complete coverage, so it is carried straight through the monthly means
# and deliberately kept out of the PCA imputation below.
######


rm(list = ls())
set.seed(202506) # for imputaion later
source('./R/utils.R')
library(dplyr)
library(lubridate)

# data are available to download from the raw csv option 
# on the CDMO website
# this script will read and process

path = '~/Library/CloudStorage/Box-Box/TGCRC Plankton Food Webs/Data/NERRFC/SWMP_TRIP_RECS'

all_files <- dir(path, full.names = TRUE)
years <- c(2011:2021)
raw_files <- list()
for(i in 1:length(years)) {
  file <- grep(years[i], all_files, value = TRUE)
  raw_files[[as.character(years[i])]] <- read.csv(file)

  if(any(grepl("TSS",names(raw_files[[as.character(years[i])]])))) {
    raw_files[[as.character(years[i])]] <- raw_files[[as.character(years[i])]] |> 
      select(-c(TSS, F_TSS))
  }
}


raw_df <- do.call(rbind, raw_files)

# # CHECK WITH SOP for QAQC DATA
# raw_df$F_PO4F |> table()
# raw_df$F_NH4F |> table()
# raw_df$F_NO23F |> table()
# raw_df$F_SALT_N |> table()
# raw_df$F_WTEM_N |> table()

# set below range to 0
raw_df$PO4F[grepl("<-4>",raw_df$F_PO4F)] <- 0
raw_df$NH4F[grepl("<-4>",raw_df$F_NH4F)] <- 0
raw_df$NO23F[grepl("<-4>",raw_df$F_NO23F)] <- 0

raw_df$SiO4F[grepl("<-4>", raw_df$F_SiO4F)] <- 0

raw_df$SiO4F[grepl("<-3>", raw_df$F_SiO4F)] <- 0
raw_df$PO4F[grepl("<-3>", raw_df$F_PO4F)] <- NA
raw_df$NH4F[grepl("<-3>", raw_df$F_NH4F)] <- NA
raw_df$NO23F[grepl("<-3>", raw_df$F_NO23F)] <- NA
raw_df$SiO4F[grepl("<-3>", raw_df$F_SiO4F)] <- NA


raw_df$Date <- raw_df$DateTimeStamp |> 
  as.POSIXct(format = '%m/%d/%Y %H:%M') |> 
  as.Date()

# remove funky nitrogen observation
raw_df$NO23F[which(
  grepl('sc',raw_df$Station.Code) & raw_df$Date %in% as.Date(c("2016-05-05",'2016-05-06'))
)] <- NA

raw_df$yearmo <- make_yearmo(raw_df, 'Date')

raw_df$sampling_site <- substr(raw_df$Station.Code, 4,5) |> toupper()
mo_nut <- raw_df |> 
  select(sampling_site, yearmo, P = PO4F, NH4 = NH4F, N = NO23F, SiOH = SiO4F, CHLA_N, temp = WTEM_N, sal = SALT_N) |>
  group_by(sampling_site, yearmo) |> 
  summarize(
    across(everything(), \(x) mean(x, na.rm = TRUE))
  ) |> 
  ungroup()

# for(col in names(mo_nut)) {print(col);print(any(is.na(mo_nut[[col]])))}
# there are some missing values for sal and temp

library(missMDA)

fix_df <- cbind(year(mo_nut$yearmo), month(mo_nut$yearmo), mo_nut |> 
  select(
    temp, sal, P, NH4, N, CHLA_N, SiOH
  )
)

fix_pca <- imputePCA(fix_df)

# plot(fix_pca$completeObs[,3], fix_pca$fittedX[,3], col = is.na(mo_nut$temp)+1)
# plot(fix_pca$completeObs[,4], fix_pca$fittedX[,4], col = is.na(mo_nut$sal)+1)


mo_nut$temp[is.na(mo_nut$temp)] <- fix_pca$completeObs[,3][is.na(mo_nut$temp)]
mo_nut$sal[is.na(mo_nut$sal)] <- fix_pca$completeObs[,4][is.na(mo_nut$sal)]

saveRDS(
  mo_nut,
  './data/01a-swmp_wq_data.rds'
)
