rm(list= ls())
library(rstan)
library(dplyr)
library(tidyr)
library(ggplot2)

source('./R/utils-mod_making.R')


etx <- readRDS("./data/02-full_merged.RDS")

taxa_group <- 'diatom'
niter = 8000
nburn = 1000
nchains = 3


core_data <- prep_component_data(taxa_group)

###################
# MARK: MODEL CONSTRUCTION 
###################

# region \- make predictors ------------
counts_poss_preds <- c('sal','t','do','turb','P','NH4','N','windspeed')
indv_poss_preds <- c('sal','t','P','NH4','N')


# region \- possible preds --------------------------------

poss_full <- list(
  full = list(
    count_preds = counts_poss_preds,
    indv_preds = indv_poss_preds
  ),
  alt = list(
    count_preds = c('sal','t','windspeed','P','NH4','N'),
    indv_preds  = indv_poss_preds
  ),
  env_only = list(
    count_preds = c('sal','t','do','turb', 'windspeed'),
    indv_preds = c('sal','t','P','NH4','N')
  ),
  red1 = list(
    count_preds = c('sal','t','turb', 'windspeed'),
    indv_preds = c('sal','t','P','NH4','N')
  ),
  red2 = list(
    count_preds = c('sal','t','windspeed'),
    indv_preds = c('sal','t','P','NH4','N')
  ),
  red3 = list(
    count_preds = c('sal','windspeed'),
    indv_preds = c('sal','t','P','NH4','N')
  ),
  red4 = list(
    count_preds = c('sal'),
    indv_preds = c('sal','t','P','NH4','N')
  )
)
poss_matrix <- list()
for(name in names(poss_full)) {
  poss_matrix[[name]] <- prep_pred_matrix(
    poss_full[[name]]$count_preds,
    poss_full[[name]]$indv_preds,
    core_data
  )
}

no_indv_preds <- lapply(
  names(poss_matrix),
  function(x) 
    return(list(
      X_counts = poss_matrix[[x]]$X_counts,
      X_indv = matrix(1, nrow(core_data$indv_mes), 1),
      pred_indv = matrix(1, nrow(core_data$all_conc), 1) # take this matrix as indv if noindv is selected
    ))
)
names(no_indv_preds) <- paste0(names(poss_matrix), '_noindv')

# region \-\- all together -------------
all_preds <- c(poss_matrix, no_indv_preds)

 
main_model <- stan_model(file = './stan/gen-hierachical_mixture_model.stan')

mod_list <- list()
pb = txtProgressBar(min = 0, max = length(names(all_preds)))
for(name in names(all_preds)) {

  progress = which(names(all_preds) == name)
  mod_list[[name]] <- fit_component_mod(
    all_preds[[name]],
    core_data,
    main_model,
    chains = nchains,
    iter = niter,
    warmup = nburn,
    cores = 7
  )
  cat(
    "\n\n\n",'---------------------------','\n',
    "Finished Model:", name, 'at', as.character(Sys.time()),
    '\n', "---------------------------",'\n\n\n'
  )
  setTxtProgressBar(pb, progress)
}


########
# MARK: Diagnositcs
########


all_rmspe <- sapply(
  names(mod_list),
  function(x) mod_list[[x]] |> rstan::extract('RMSPE')
)

rmspe_summary <- sapply(all_rmspe, quantile, probs = c(0.025,0.5,0.975))
write.csv(rmspe_summary,paste0('./output/03-',taxa_group,'-rmspe.csv'))

plot_df <- sapply(
  all_rmspe,
  function(x) x[which(x < quantile(x, 0.975) & x > quantile(x, 0.025))]
) |> 
  as.data.frame()

min_name <- all_rmspe |> 
  sapply(mean) |> 
  which.min() |> 
  names()

ggplot(plot_df |> 
  pivot_longer(
    cols = everything(),names_to = 'Model',values_to = 'RMSPE'
  )) +
  geom_violin(
    aes(
      x = Model,
      y = RMSPE
    )
  ) +
  geom_point(
    aes(
      x = Model,
      y = mean(RMSPE)
    )
  )+
  labs(subtitle = paste0('Best Model is ',min_name))+
  theme_classic()+
  theme(axis.text = element_text(angle = 45))


ggsave(paste0('./output/03-best_fit-', taxa_group,'.png'))