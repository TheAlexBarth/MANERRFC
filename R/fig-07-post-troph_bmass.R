rm(list = ls())
source('./R/utils.R')
post <- readRDS("./data/04-post_micro.RDS")

conc_raw <- readRDS("./data/02-full_merged.rds")$conc

post$alpha |> dim()


# region CALCULATE Lambda exp


lambda_list <- list()
for(troph in troph_factors) {
  lambda_list[[troph]] <- list()
  troph_idx <- which(levels(troph_factors) == troph)
  for(site in site_factors) {
    site_idx <- which(levels(site_factors) == site)

    data_idx <- which(conc_raw$functional_role == troph & conc_raw$sample_site == site)
    lambda_list[[troph]][[site]] <- as.matrix(post$preds[data_idx,]) %*% t(post$alpha[,troph_idx, site_idx, ]) |> 
      exp() |> 
      t() |> 
      summarize_pred()

    

    lambda_list[[troph]][[site]]$date <- conc_raw$yearmo[data_idx]

  }
  lambda_list[[troph]] <- EcotaxaTools::list_to_tib(lambda_list[[troph]], 'site')
}

# endregion