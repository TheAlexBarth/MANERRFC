######
# Silicate (SiOH) Formatting
#
# Added for revision per reviewer request. SiO4F has near-complete monthly
# coverage at every site across the full plankton sampling window (see
# coverage check below), so no imputation is needed here (unlike pipe-01a's
# PCA imputation for temp/sal).
######

rm(list = ls())
source('./R/utils.R')
library(dplyr)
library(lubridate)

# data are available to download from the raw csv option
# on the CDMO website
# will need local adjustment

path <- '~/Library/CloudStorage/Box-Box/TGCRC Plankton Food Webs/Data/NERRFC/SWMP_TRIP_RECS'

all_files <- dir(path, full.names = TRUE)
years <- c(2011:2021)
raw_files <- list()
for(i in 1:length(years)) {
  file <- grep(years[i], all_files, value = TRUE)
  raw_files[[as.character(years[i])]] <- read.csv(file) |>
    select(Station.Code, DateTimeStamp, SiO4F, F_SiO4F, PO4F, F_PO4F, NH4F, F_NH4F, NO23F, F_NO23F)
}

raw_df <- do.call(rbind, raw_files)

# set below range to 0 / rejected to NA, same convention as pipe-01a
raw_df$SiO4F[grepl("<-4>", raw_df$F_SiO4F)] <- 0
raw_df$SiO4F[grepl("<-3>", raw_df$F_SiO4F)] <- NA
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
raw_df$sampling_site <- substr(raw_df$Station.Code, 4, 5) |> toupper()

#confirm that there aren't too many NAs
raw_df |> apply(1, function(x) {
    any(is.na(x[c('SiO4F','PO4F','NH4F','NO23F')]))
  }) |>
  sum() / nrow(raw_df)
# it's less than 5% ok cool

mo_si <- raw_df |>
  select(sampling_site, yearmo, SiOH = SiO4F) |>
  group_by(sampling_site, yearmo) |>
  summarize(
    SiOH = mean(SiOH, na.rm = TRUE),
    .groups = 'drop'
  )

saveRDS(
  mo_si,
  './data/01e-silicate.rds'
)

# region \- coverage check against plankton sampling occasions -----------------
# confirms SiOH coverage is high enough at the site x yearmo combinations
# already used in the plankton dataset before this predictor is joined in
# downstream (pipe-02-data_merger.R, 03-/04- regression scripts)
etx <- readRDS('./data/00-ecotaxa_full.rds')
plankton_occasions <- etx$conc |>
  distinct(sample_site, yearmo) |>
  rename(sampling_site = sample_site) |>
  mutate(sampling_site = as.character(sampling_site))

coverage <- plankton_occasions |>
  left_join(mo_si, by = c('sampling_site', 'yearmo'))

cat('Plankton site-months:', nrow(plankton_occasions), '\n')
cat('...with SiOH available:', sum(!is.na(coverage$SiOH)), '\n')
cat('Overall coverage:', round(100 * mean(!is.na(coverage$SiOH)), 1), '%\n')
# endregion -----------------------
