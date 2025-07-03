rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

source("./R/utils.r")

# region Prep data -----------------------
post <- readRDS('./data/04-post_micro.RDS')
raw <- readRDS("./data/02-full_merged.rds")

marg_list <- list()


# region Effect table -----------------------

group_levels <- list(
  'troph' = levels(troph_factors),
  'site' = site_factors,
  'w' = names(post$preds[,-1])
)

alpha_df <- post_arr_to_df(post$alpha[,,,-1], c('troph','site','w'), dim_levels = group_levels)


alpha_sum <- alpha_df |> 
  group_by(troph, site, w) |> 
  summarize(
    low = quantile(est, probs = c(0.025)),
    mid = median(est),
    high = quantile(est, probs = c(0.975)),
    pos_prob = mean(est > 0)
  )



effect_table <- alpha_sum |> 
  select(w, mid, low, high, pos_prob) |> 
  mutate(
    `Post. Mean` = sprintf('%.2f', mid),
    `95% Cred.Intv` = paste0("(", sprintf("%.2f", low), ", ", sprintf("%.2f", high), ")"),
    `Pr(>0)` = sprintf("%.2f", pos_prob)
  ) |> 
  select(`Trophic Role` = troph, `Site` = site, Variable = w, `Post. Mean`, `95% Cred.Intv`, `Pr(>0)`)

  write.csv(effect_table, './output/s-troph_post_summary.csv', row.names = FALSE)
  

# endregion -----------------------

# region seasonal
# for(troph in troph_factors) {
#   marg_list[['month']][[troph]] <- list()
#   troph_idx <- which(levels(troph_factors) == troph)
#   mean_bmass <- exp(post$eta_r[,troph_idx] + (post$sigma_r[,troph_idx]^2/2)) # pgC
#   for(site in site_factors) {
#     site_idx <- which(levels(site_factors) == site)

#     sin_alpha <- post$alpha[,troph_idx, site_idx, 2]
#     cos_alpha <- post$alpha[,troph_idx, site_idx, 3]

#     amp <- sqrt(sin_alpha^2 + cos_alpha^2)
#     phase <- atan2(sin_alpha, cos_alpha)

#     lambda_pred <- sapply(1:12, function(t){
#       post$alpha[,troph_idx, site_idx, 1] +
#         amp* sin(2*pi*t/12 + phase)
#     })

#     pconc <- t(exp(lambda_pred)*1e6) * rep(mean_bmass/1e9, each = ncol(lambda_pred))

#     marg_list[['month']][[troph]][[site]] <- pconc |> 
#       t() |> 
#       summarize_pred()
#     marg_list[['month']][[troph]][[site]]$month <- 1:12
#   }
#   marg_list[['month']][[troph]] <-   marg_list[['month']][[troph]] |> 
#     EcotaxaTools::list_to_tib('site')
# }

# marg_list[['month']] <- marg_list[['month']] |> 
#   EcotaxaTools::list_to_tib('troph')

# region All Vars -----------------------
for(var in names(post$preds[,-1])) {
  marg_list[[var]] <- list()
  var_idx <- which(names(post$preds) == var)
  for(troph in troph_factors) {
    if(var %in% size_factors & grepl('auto', troph)) {
      next
    }
    marg_list[[var]][[troph]] <- list()
    troph_idx <- which(levels(troph_factors) == troph)

    for(site in site_factors) {
      site_idx <- which(levels(site_factors) == site)
      
      sim_range <- seq(
        min(raw$conc[[var]]),
        max(raw$conc[[var]]),
        length.out = 100
      )

      if(var %in% c('N',"NH4","P")) {
        sim_scale <- scale(
          log(sim_range +1e-5), 
          center = mean(log(raw$conc[[var]]+1e-5)), 
          scale = sd(log(raw$conc[[var]]+1e-5)))
      } else {
        sim_scale <- scale(
          sim_range,
          center = mean(raw$conc[[var]]),
          scale = sd(raw$conc[[var]])
        )
      }

      lambda_pred <- post$alpha[,troph_idx, site_idx, 1] +
          as.matrix(sim_scale) %*% post$alpha[, troph_idx, site_idx, var_idx]

      pconc <- t(exp(lambda_pred)*1e6) * rep(mean_bmass/1e9, each = nrow(lambda_pred))

      marg_list[[var]][[troph]][[site]] <- (exp(lambda_pred)*1e6) |> 
        t() |> 
        summarize_pred()
      marg_list[[var]][[troph]][[site]][[var]] <- sim_range
    }
    marg_list[[var]][[troph]] <-   marg_list[[var]][[troph]] |> 
      EcotaxaTools::list_to_tib('site')
  }
  marg_list[[var]] <- marg_list[[var]] |> 
    EcotaxaTools::list_to_tib('troph')
}



# endregion -----------------------


# region Plot -----------------------

marg_plotter <- function(var) {
    p = ggplot(
      marg_list[[var]] |> 
        mutate(
          troph = factor(troph, levels = troph_factors)
        )
    ) +
      geom_ribbon(
        aes(
          x = .data[[var]],
          ymin = log(low.95),
          ymax = log(high.95),
          fill = troph
        ),
        alpha = 0.5
      ) +
      # geom_ribbon(
      #   aes(
      #     x =  marg_list[[var]][[var]],
      #     ymin = low.75,
      #     ymax = high.75,
      #     fill = troph
      #   ),
      #   alpha = 0.5
      # ) +
      # geom_ribbon(
      #   aes(
      #     x =  marg_list[[var]][[var]],
      #     ymin = low.50,
      #     ymax = high.50,
      #     fill = troph
      #   ),
      #   alpha = 0.5
      # ) +
    
      facet_grid(
        site ~ troph
      ) +
      scale_fill_manual(
        values = troph_cols
      )+
      guides(fill = 'none')+
      labs(x = var, y = 'Log(Num per L)')+
      theme_pubclean() +
      theme(strip.background = element_rect(fill = 'transparent'))
  return(p)
}


pdf('./output/s-post_troph_marg.pdf')
for(var in names(marg_list)) {
  print(marg_plotter(var))
}
dev.off()
# endregion -----------------------