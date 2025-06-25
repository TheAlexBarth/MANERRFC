######
# Environmental Formatting
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
years <- c(2014:2021)
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

raw_df$PO4F[grepl("<-3>", raw_df$F_PO4F)] <- NA
raw_df$NH4F[grepl("<-3>", raw_df$F_NH4F)] <- NA
raw_df$NO23F[grepl("<-3>", raw_df$F_NO23F)] <- NA


raw_df$Date <- raw_df$DateTimeStamp |> 
  as.POSIXct(format = '%m/%d/%Y %H:%M') |> 
  as.Date()
raw_df$yearmo <- make_yearmo(raw_df, 'Date')

raw_df$sampling_site <- substr(raw_df$Station.Code, 4,5) |> toupper()
mo_nut <- raw_df |> 
  select(sampling_site, yearmo, P = PO4F, NH4 = NH4F, N = NO23F, CHLA_N, temp = WTEM_N, sal = SALT_N) |> 
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
    temp, sal, P, NH4, N, CHLA_N
  )
)

fix_pca <- imputePCA(fix_df)

mo_nut$temp[is.na(mo_nut$temp)] <- fix_pca$completeObs[,3][is.na(mo_nut$temp)]
mo_nut$sal[is.na(mo_nut$sal)] <- fix_pca$completeObs[,4][is.na(mo_nut$sal)]

saveRDS(
  mo_nut,
  './data/01a-swmp_wq_data'
)
