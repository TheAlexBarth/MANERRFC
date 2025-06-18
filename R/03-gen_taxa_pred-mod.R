rm(list= ls())
library(cmdstanr)
library(dplyr)
library(tidyr)
library(ggplot2)

source('./R/utils-mod_making.R')


etx <- readRDS("./data/02-full_merged.RDS")

taxa_group <- 'mz'
log_file <- './logs/20250617-mz_run.log'
file.create(log_file)




###################
# MARK: MODEL CONSTRUCTION 
###################



# region \-\- all together -------------
 
main_model <- cmdstan_model('./stan/gen-hierachical_mixture_model.stan')

bmass_priors <- list(
  b_mu = 8,
  b_sig = 3,
  tau_mu = 1,
  tau_sig = 1.5
)


# region \- make predictors ------------
counts_poss_preds <- c('sal','t','do','turb','CHLA_N','windspeed', "TotalPAR", 'P','NH4','N')

cat(
  '\n','Launching Model Selection','\n',
  sep = "\n", file = log_file,
  append = TRUE
)
full_data <- prep_component_data(taxa_group, counts_poss_preds)
full_mod <- fit_mod(
  data_list = full_data,
  count_preds = counts_poss_preds,
  bmass_priors,
  main_model,
  iter = 7000,
  chains = 10,
  parallel_chains = 10,
  iter_warmup = 2000,
  thin = 5,
  refresh = 0
)
# make sure to check before running that the full model passes post-pred check
cur_pred = counts_poss_preds
cur_mod = full_mod
cur_rmspe <- cur_mod$draws('RMSPE',,'df')$RMSPE
improve = TRUE

while(improve & length(cur_pred)>1) {
  improve <- FALSE
  cat(
    "",paste0('Current Pars: ', paste0(cur_pred, collapse = ', ')),
    sep = "\n", file = log_file,
    append = TRUE
  )
  for(i in rev(seq_along(cur_pred))) {
    drop_pred <- cur_pred[i]

    cat(
      "",'=================',paste0('Checking: ', drop_pred), '=================', 
      sep = "\n", file = log_file,
      append = TRUE
    )

    can_data <- prep_component_data(taxa_group, cur_pred[-i])
    can_mod <- fit_mod(
      data_list = can_data,
      count_preds = cur_pred[-i],
      bmass_priors,
      main_model,
      iter = 4000,
      chains = 5,
      parallel_chains = 5,
      iter_warmup = 2000,
      thin = 2,
      refresh = 0
    )
    # region \-\- check bayes p value ----------
    pdf <- can_mod$draws(c('MSE_obs','MSE_til'), format = 'df')
    p <- mean(pdf$MSE_obs > pdf$MSE_til)
    if(p > 0.95 | p < 0.05) {
      p_accept <- FALSE
    } else {
      p_accept <- TRUE
    }
    # endregion

    # region \-\- check RMSE ------------
    can_rmspe <- can_mod$draws('RMSPE',,'df')$RMSPE
    p_better <- mean(can_rmspe <= cur_rmspe)
    # endregion

    # region \-\- check effect -------------
    # note this checks the currrent model not candidate!
    betas <- cur_mod$draws('beta_count', format = 'df')
    no_effect <- betas |> 
      select(grep(paste0(i+1,"]$"),names(betas))) |> 
      apply(2, quantile, probs = c(0.025,0.975)) |> 
      overlaps_zero() |> 
      all()
    # endregion

    if(p_better > 0.40 & p_accept & no_effect) {
      improve <- TRUE
      cur_mod <- can_mod
      cur_pred <- cur_pred[-i]
      cur_rmspe <- can_rmspe
      cat(
        '====================',
        paste0('Dropping: ', drop_pred), '====================',
        sep = "\n", file = log_file,
        append = TRUE
      )
      break
    } else {
      # log message:
      if(p_better <= 0.40) {
        cat(
          '====================',
          paste0('Keeping: ', drop_pred, ' Worsened RMPSE beyond threshold'), '====================',
          sep = "\n", file = log_file,
          append = TRUE
        )
      }
      if(!p_accept) {
        cat(
          '====================',
          paste0('Keeping: ', drop_pred, ' PPC fail'), '====================',
          sep = "\n", file = log_file,
          append = TRUE
        )
      }
      if(!no_effect) {
        cat(
          '====================',
          paste0('Keeping: ', drop_pred, ' Effect at Some Level'), '====================',
          sep = "\n", file = log_file,
          append = TRUE
        )
      }
    }
  }
}
cat(
  '====================',
  paste0('Model Finalized: ', paste0(cur_pred, collapse = ', ')),
  '====================',
  sep = "\n", file = log_file,
  append = TRUE
)

# manual tests:

modterms_1 = c("sal", "t", "do", "P", "NH4")

modterms_2 = c("sal", "t", "do", "P", "NH4")


# # consider a manual check for dropping all nutrient parameters
# #
# # manual check for no nutrient params:
# can_preds <- c('sal','t','do','turb','windspeed')
# can_data <- prep_component_data(taxa_group, can_preds)
# can_mod <- fit_mod(
#   data_list = can_data,
#   count_preds = can_preds,
#   bmass_priors,
#   main_model,
#   iter = 4000,
#   chains = 5,
#   parallel_chains = 5,
#   iter_warmup = 2000,
#   thin = 2,
#   refresh = 0
# )
# did not improve results...

ggplot() +
  geom_density(
    aes(x = can_rmspe),
    fill = 'blue', alpha = 0.25
  ) +
  geom_density(
    aes(x = cur_rmspe),
    fill = 'red', alpha = 0.25
  )
