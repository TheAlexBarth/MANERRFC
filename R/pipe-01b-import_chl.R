# rscript to format chlorophyll data from CDMO
# need to get messy excel files from CDMO then this can read them

rm(list = ls())

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)

path = '~/Library/CloudStorage/Box-Box/TGCRC Plankton Food Webs/Data/NERRFC/size_frac_chl'

all_files <- dir(path, full.names = TRUE)


extract_date <- function(x) {
  x <- as.character(x)
  if(!is.na(as.numeric(x))) {
    return(as.Date(as.numeric(x), origin = '1899-12-30'))
  } else {
    x <- str_extract(x, "\\d{1,2}/\\d{1,2}/\\d{2,4}")
    if (!is.na(x)){
      return(as.Date(x, format = "%m/%d/%Y"))
    } else {
      return(NA)
    }
  }
}

full_list <- list()
for(i in 1:length(all_files)) {
  year = c(2014:2021)[i]
  file = all_files[i]
  raw <- read_excel(file, sheet = 'Means', col_names = FALSE, col_types = "text")

  res <- list()
  cur_sites <- NULL
  cur_date <- NA
  row_num <- 1

  while(row_num <= nrow(raw)) {
    row <- raw[row_num,]

    # get date if possible
    cell_one <- as.character(row[[1]])

    if(is.na(cell_one)) {
      row_num = row_num + 1
      next
    }

    date_try <- extract_date(cell_one) |> suppressWarnings()
    if(!is.na(date_try)) {
      cur_date <- date_try
      row_num <- row_num + 1
      next
    }

    # if collected row exists, get new dates
    if(toupper(cell_one) == "COLLECTED") {
      collect_dates <- raw[row_num, -1] |> 
        sapply(extract_date) |> 
        as.Date()
      cur_date <- collect_dates
      row_num <- row_num + 1
      next
    }

    # measurement block

    if(cell_one %in% c("GF/F", "5-20 um", "<20.0 um")) {
      type = cell_one
      mean_row <- row_num + 2

      means <- as.numeric(raw[mean_row,-1])

      res[[length(res) + 1]] <- data.frame(
        year = year,
        date = cur_date,
        site = c('AB','CE','CW','MB','SC'),
        meas = type,
        value = means
      )
    }
    row_num = row_num + 1

  }

  full_list[[i]] <- do.call(what = rbind, res)
}

full_data <- do.call(what = rbind, full_list)

full_data$meas <- full_data$meas |> 
  sapply(function(x) {
    switch(x,
      'GF/F' = 'Whole',
      '5-20 um' = 'five',
      '<20.0 um' = 'twenty'
    )
  })


wide_data <- full_data |> 
  pivot_wider(values_from = value, names_from = meas)

# conceptually these cannot be less than 0
pos_only <- function(val) {
  ifelse(val > 0, val, 0)
}


wide_data$micro <- pos_only(wide_data$Whole - wide_data$twenty) # subtract under twenty
wide_data$nano <- wide_data$five # this is nano section
wide_data$pico <- pos_only(wide_data$twenty - wide_data$five) # under twenty minus 5-20

good_data <- wide_data |> 
  select(
    year, date,
    site, micro, nano, pico
  )

good_data$yearmo <- as.Date(paste(year(good_data$date), month(good_data$date), '01',sep = '-'))

saveRDS(good_data,'./data/01b-size_frac_chl.RDS')
