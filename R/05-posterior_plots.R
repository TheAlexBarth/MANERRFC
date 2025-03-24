rm(list = ls())
library(ggplot2)
library(tidyr)
library(dplyr)



#############
# MARK: CHOOSE TAXA 
#############

taxa = 'diatom'


##############################
# MARK: LOAD DATA
##############################


full <- readRDS('./data/02-full_merged.rds')
mod = readRDS(paste0('./data/04-mod_output-',taxa,'.RDS'))

source('./R/utils.R')


####################################
# MARK: RUN TO SET UP 
####################################



# region \- background set up ----------------------------

# this is slow so best to do once outside the function
betas <- rstan::extract(mod$mod, 'beta_count')$beta_count

beta_list <- list()
for(i in 1:dim(betas)[3]) {
  beta_list[[i]] <- betas[,,i]
  colnames(beta_list[[i]]) <- mod$group_names
}
if(length(beta_list) == length(mod$count_preds)+1) {
  names(beta_list) <- c('int',mod$count_preds)
} else {
  stop("Beta names are wrong")
}

make_effect_plot <- function(param) {
  if(!(param) %in% names(beta_list)) {
    stop(paste0("Wrong param, only can use: ", paste(names(beta_list), collapse = ' or ')))
  }

  df <- beta_list[[param]] |> 
    as.data.frame() |> 
    pivot_longer(everything(), names_to = 'taxa', values_to = 'est')

  df_sum <- df |> 
    group_by(taxa) |> 
    summarize(
      mean = mean(est),
      lower = quantile(est, probs = 0.025),
      upper = quantile(est, probs = 0.975)
    )
  
  p1 = ggplot() +
    geom_violin(
      data = df,
      aes(
        x = taxa,
        y = est,
        fill = taxa
      )
    ) +
    geom_hline(yintercept = 0, col = 'red')+
    geom_point(
      data = df_sum,
      aes(
        x = taxa,
        y = mean
      ),
      col = 'black',
      size = 3
    ) +
    geom_segment(
      data = df_sum,
      aes(
        x = taxa,
        y = lower,
        yend = upper
      ),
      inherit.aes = FALSE,
      linewidth = 0.5,
      col = 'black'
    )+
    theme_classic() +
    labs(x = "",fill = "")
    theme(axis.text.x = element_blank())
  
    return(p1)
}

####################
# MARK: MAKE PLOTS HERE 
###################


make_effect_plot('sal') + theme(axis.text.x = element_text(angle = 45))



sal_plot <- make_effect_plot('sal')

sal_plot +
  labs(subtitle = 'Salinity')


t_plot <- make_effect_plot('t')

wind_plot <- make_effect_plot('windspeed')

