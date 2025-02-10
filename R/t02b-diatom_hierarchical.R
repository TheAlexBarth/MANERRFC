rm(list = ls())
library(EcotaxaTools) # this is my package
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(rstan)
library(lubridate)
source('./R/utils.R')

conc <- readRDS('./data/t01-copano_concentrations.rds')
name_list = readRDS('./data/t01-taxa_names.rds')
indv <- readRDS('./data/t01-copano_individals.rds')
set.seed(1200)

# trim to diatoms
conc <- conc |> 
  filter(taxa %in% name_list$diatom)

# flip taxa names 
conc$group <- conc$taxa |> 
  sapply(
    function(x) names(name_list$diatom)[which(name_list$diatom == x)]
  ) |> 
  as.factor()


#########################
# MARK: Preliminary plots
#########################
# g_sal <- list()
# for(g in unique(conc$group)) {
#   print(g)
#   g_sal[[g]] <- ggplot() +
#     geom_point(
#       data = conc |> 
#         filter(group == g),
#       aes(
#         x = sal,
#         y = um3
#       )
#     ) +
#     labs(subtitle = g)
# }
# for(i in 1:length(g_sal)) {
#   quartz()
#   print(g_sal[[i]])
# }


#####
# MARK: data prep 
######

all_diat <- conc |>
  group_by(sample_id, group) |> 
  summarize(count = sum(count), num_L = sum(num_L), um3_L = sum(um3_L)) |> 
  left_join(
    conc |> 
    select(-c(count, num_L, um3,um3_L, taxa, group)) |> 
      unique(),
    by = "sample_id"
  )

# impute to mean
all_diat$sal[which(is.na(all_diat$sal))] <- mean(all_diat$sal, na.rm = T)
all_diat$t[which(is.na(all_diat$t))] <- mean(all_diat$t, na.rm = T)

indv_metrics <- indv |> 
  filter(taxo_name %in% name_list$diatom)

indv_metrics$group <- indv_metrics$taxo_name |> 
  sapply(
    function(x) names(name_list$diatom)[which(name_list$diatom == x)]
  ) |> 
  as.factor()



####
# MARK: MODEL 
####


# region \- stan fit ----------
indv_metrics$sal[is.na(indv_metrics$sal)] <- mean(indv_metrics$sal, na.rm = T)

# don't use mix-multiple since there is no effect of temperature
pois_data <- list(
  N_obs = nrow(all_diat),
  N_mes = nrow(indv_metrics),
  N_groups = length(unique(all_diat$group)),
  n = all_diat$count,
  group = as.numeric(all_diat$group),
  img_vol = all_diat$img_vol/1000, # convert to L,
  sal = all_diat$sal,
  group_counts = as.numeric(all_diat$group),
  group_bio = as.numeric(indv_metrics$group),
  log_b = log(indv_metrics$um3),
  sal_b = indv_metrics$sal,
  temp_b = indv_metrics$t,
  obs_biomass_conc = all_diat$um3_L
)

pois_fit <- stan(
  file = './stan/03b-mix-sal_only-hierachical.stan',
  data = pois_data,
  chains = 4, iter = 3000, warmup = 500, cores = 14
)


####################
# MARK: POSTEIOR EVALS
####################


# region \- effect level -------------------------------
sal_pred = seq(min(all_diat$sal, na.rm = T), max(all_diat$sal, na.rm = T), by = 0.1)


# predict lambda
g_pred_list <- list()
lambda_fit <- list()
n_pred <- list()
for(l in 1:length(levels(all_diat$group))) {
  lambda_fit[[l]] <- post_predict(
    extract(pois_fit),
    list(
      sal = sal_pred
    ),
    paste0(
      'exp(beta0[,',l,'] + beta1[,',l,'] * sal + log(0.1129333/1000))'
    ),
    n_draws = 2000
  )

  n_pred[[l]] <- matrix(
    rpois(length(lambda_fit[[l]]), lambda = lambda_fit[[l]]),
    nrow = nrow(lambda_fit[[l]]),
    ncol = ncol(lambda_fit[[l]])
  )
}


# Predict Biomass
err_var <- extract(pois_fit, pars = 'v')$v
eta_pred <- list()
d_pred <- list()
for(l in 1:length(levels(all_diat$group))) {
  eta_pred[[l]] <- post_predict(
    extract(pois_fit),
    list(
      sal = sal_pred
    ),
    paste0(
      'a0[,',l,'] + a1[,',l,'] * sal'
    ),
    n_draws = 2000
  )
 
  d_pred[[l]] <- apply(
    eta_pred[[l]], 2, 
    function(x) rnorm(nrow(eta_pred[[l]]), mean = x, sd = err_var)
  ) |> 
    exp()

  g_pred_list[[l]] = d_pred[[l]] * (n_pred[[l]]/mean(all_diat$img_vol/1000))
}



# region \- TS PRED ---------------------

# super messy need to fix, add code

wq = wq_daysum |> 
    filter(sample_site %in% c('CE', 'CW')) |> 
    group_by(Date) |> 
    summarise(sal = mean(Sal))

wq = wq[!is.na(wq$sal),]


# predict lambda
g_ts_list <- list()
lamba_ts <- list()
n_ts <- list()
for(l in 1:length(levels(all_diat$group))) {
  lamba_ts[[l]] <- post_predict(
    extract(pois_fit),
    list(
      sal = wq$sal
    ),
    paste0(
      'exp(beta0[,',l,'] + beta1[,',l,'] * sal + log(0.1129333/1000))'
    ),
    n_draws = 2000
  )

  n_ts[[l]] <- matrix(
    rpois(length(lamba_ts[[l]]), lambda = lamba_ts[[l]]),
    nrow = nrow(lamba_ts[[l]]),
    ncol = ncol(lamba_ts[[l]])
  )
}


# Predict Biomass
err_ts <- extract(pois_fit, pars = 'v')$v
eta_ts <- list()
d_ts <- list()
for(l in 1:length(levels(all_diat$group))) {
  eta_ts[[l]] <- post_predict(
    extract(pois_fit),
    list(
      sal = wq$sal
    ),
    paste0(
      'a0[,',l,'] + a1[,',l,'] * sal'
    ),
    n_draws = 2000
  )
 
  d_ts[[l]] <- apply(
    eta_ts[[l]], 2, 
    function(x) rnorm(nrow(eta_ts[[l]]), mean = x, sd = err_ts)
  ) |> 
    exp()

  g_ts_list[[l]] = d_ts[[l]] * (n_ts[[l]]/mean(all_diat$img_vol/1000))
}



############################
# MARK: GROUP PLOTTING ###########
############################
 
solo_pennate <- g_pred_list[[7]] |> 
  summarize_pred()

# region \- pennate conc - - - - - - - - - - - - 
ggplot() +
  geom_point(
    data = all_diat |> 
      filter(group == 'solo_pennate'),
    aes(
      x = sal,
      y = um3_L,
      color = gg_cbb_col(7)[7]
    )
  ) +
  geom_line(
    aes(
      x = sal_pred,
      y = solo_pennate$mean
    )
  ) +
  geom_ribbon(
    aes(
      x = sal_pred,
      ymin = solo_pennate$low,
      ymax = solo_pennate$high
    ),
    alpha = 0.25
  ) +
  theme_pubclean()+
  theme(legend.position = 'none')



# region \- pennate biomass  - - - - - - 
pen_bmass <- d_pred[[7]] |> 
  summarize_pred() |> 
  mutate(sal = sal_pred)

ggplot() +
  geom_point(
    data = indv_metrics |> 
      filter(group == "solo_pennate"),
    aes(
      x = sal,
      y =log(um3)
    )
  ) + 
  geom_line(
    data = pen_bmass,
    aes(
      x = sal,
      y = log(mean)
    )
  ) +
  geom_ribbon(
    data = pen_bmass,
    aes(
      x = sal,
      ymin = log(low),
      ymax = log(high)
    ),
    alpha = 0.25
  )



Thalassionema <- g_pred_list[[8]] |> 
  summarize_pred()

ggplot() +
  geom_point(
    data = all_diat |> 
      filter(group == 'Thalassionema'),
    aes(
      x = sal,
      y = um3_L
    ),
    color = gg_cbb_col(8)[8],
    size = 3
  ) +
  geom_line(
    aes(
      x = sal_pred,
      y = Thalassionema$mean
    ),
    color = gg_cbb_col(8)[8],
    linewidth = 1.5
  ) +
  geom_ribbon(
    aes(
      x = sal_pred,
      ymin = Thalassionema$low,
      ymax = Thalassionema$high
    ),
    alpha = 0.25,
    fill = gg_cbb_col(8)[8]
  ) +
  labs(y = 'Biovolume Concentration [um3/L]', x = 'Salinity')+
  theme_pubclean()+
  theme(legend.position = 'none') +
  theme(
    axis.title = element_text(face = 'bold', size = 16),
    axis.text = element_text(face = 'bold', size = 12)
  )
  


  Thalassionema <- g_pred_list[[1]] |> 
    summarize_pred()
  
  ggplot() +
    geom_point(
      data = all_diat |> 
        filter(group == 'Baxillaria'),
      aes(
        x = sal,
        y = um3_L
      ),
      color = gg_cbb_col(8)[8],
      size = 3
    ) +
    geom_line(
      aes(
        x = sal_pred,
        y = Thalassionema$mean
      ),
      color = gg_cbb_col(8)[8],
      linewidth = 1.5
    ) +
    geom_ribbon(
      aes(
        x = sal_pred,
        ymin = Thalassionema$low,
        ymax = Thalassionema$high
      ),
      alpha = 0.25,
      fill = gg_cbb_col(8)[8]
    ) +
    labs(y = 'Biovolume Concentration [um3/L]', x = 'Salinity')+
    theme_pubclean()+
    theme(legend.position = 'none') +
    theme(
      axis.title = element_text(face = 'bold', size = 16),
      axis.text = element_text(face = 'bold', size = 12)
    )
    
  


# region \- mean trend
# mean_trend <- post_predict(
#   extract(pois_fit),
#   list(
#     sal = sal_pred
#   ),
#   paste0(
#     'exp(mu0+ mu1 * sal) * exp(theta0 + theta1 * sal)'
#   ),
#   n_draws = 2000
# ) |> 
#   summarize_pred()

# ggplot() +
#   geom_point(
#     data = all_diat,
#     aes(
#       x = sal,
#       y = um3_L,
#       color = group
#     )
#   ) +
#   geom_line(
#     aes(
#       x = sal_pred,
#       y = mean_trend$mean
#     )
#   ) +
#   geom_ribbon(
#     aes(
#       x = sal_pred,
#       ymin = mean_trend$low,
#       ymax = mean_trend$high
#     ),
#     alpha = 0.25
#   ) +
#   scale_color_manual(values = gg_cbb_col(length(levels(all_diat$group))))+
#   theme_pubclean()+
#   theme(legend.position = 'none')

############################
# MARK: TOTAL PLOTS #######
############################


# region \- summarize list fx - - - - - - - - 

summarize_list <- function(mat_list, sp_range = sal_pred){
  l = mat_list |> 
    lapply(summarize_pred) |> 
    lapply(
      function(x) mutate(x, sal = sp_range)
    ) |> 
    list_to_tib()|> 
    group_by(sal) |> 
    summarize(
      mean = sum(mean),
      low.50 = sum(low.50),
      high.50 = sum(high.50),
      low.75 = sum(low.75),
      high.75 = sum(high.75),
      low.95 = sum(low.95),
      high.95 = sum(high.95)
    ) |> 
    ungroup()
  
  return(l)
}

total_diat <- all_diat |> 
  group_by(sample_id, date) |> 
  summarize(
    biomass = sum(um3_L, na.rm = T),
    sal = mean(sal, na.rm = T)
  )


# region \- ts conc ------------------------

total_ts <- g_ts_list |> 
  summarize_list(sp_range = wq$sal)

total_plot_ts <- wq |> 
  left_join(
    total_ts,
    by = 'sal'
  )

ggplot() +
  geom_point(
    data = total_diat,
    aes(
      x = date,
      y = biomass
    ),
    size = 3, color = '#dbf6f6'
  )+
  geom_error_range(total_plot_ts$Date, total_plot_ts, '#ff9191') +
  scale_x_date(
    limits = c(min(total_plot_ts$Date), max(total_plot_ts$Date))
  )+
  labs(x = "", y = '')+
  scale_y_continuous(limits = c(0,6.5e10))+
  theme_minimal() +
  theme(
    axis.text = element_text(size = 28, color = 'white'),
    panel.grid = element_line(color = 'grey50'),
    plot.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.background = element_rect(fill = 'transparent')
  )

ggsave('./output/02b_tot-ts.pdf',
height = 6, width = 22, units = 'in')
    

# region \- total conc  - - - - - - - - - 
total_pred <- g_pred_list |> 
  summarize_list(sp_range = sal_pred)

total_diat <- all_diat |> 
  group_by(sample_id) |> 
  summarize(
    biomass = sum(um3_L, na.rm = T),
    sal = mean(sal, na.rm = T)
  )

ggplot() +
  geom_point(
    aes(
      x = total_diat$sal,
      y = total_diat$biomass
    ),
    size = 3, color = '#dbf6f6'
  )+
  geom_error_range(sal_pred, total_pred, '#ff9191')+
  labs(x = "", y = "")+
  theme_pubclean()+
  theme(
    axis.text = element_text(size = 20, color = 'white'),
    panel.grid = element_line(color = 'grey50'),
    plot.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.background = element_rect(fill = 'transparent')
  )
ggsave('./output/02b-tot_diat03-conc_pred.pdf',
height = 6, width = 6, units = 'in')



# region \- total count  - - - - - - - - - - - - -  -

total_lambda <- lambda_fit |> 
  summarize_list()

total_npred <- n_pred |> 
  summarize_list()

tobs_count <- all_diat |> 
  group_by(sample_id) |> 
  summarize(
    count = sum(count, na.rm = T),
    sal = mean(sal, na.rm = T)
  )

ggplot() +
  geom_point(
    aes(
      x = tobs_count$sal,
      y = tobs_count$count
    ),
    size = 3, color = '#dbf6f6'
  )+
  geom_error_range(sal_pred, total_npred, '#ff9191')+
  labs(x = "", y = "")+
  theme_pubclean() +
  theme(
    axis.text = element_text(size = 20, color = 'white'),
    panel.grid = element_line(color = 'grey50'),
    plot.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.background = element_rect(fill = 'transparent')
  )
ggsave('./output/02b-tot_diat01-pois_pred.pdf',
height = 6, width = 6, units = 'in')



# region \- total bmass  - - - - - - - - - - - --  -
total_d <- d_pred |> 
  lapply(summarize_pred) |> 
  lapply(
    function(x) mutate(x, sal = sal_pred)
  ) |> 
  list_to_tib()|> 
  group_by(sal) |> 
    summarize(
      mean = sum(mean),
      low.50 = sum(low.50),
      high.50 = sum(high.50),
      low.75 = sum(low.75),
      high.75 = sum(high.75),
      low.95 = sum(low.95),
      high.95 = sum(high.95)
    ) |> 
    ungroup()

total_bmass <- indv_metrics |> 
  group_by(sample_id) |> 
  summarize(
    sal = mean(sal, na.rm = T),
    um3 = mean(um3, na.rm = T)
  )


ggplot() +
  geom_point(
    aes(
      x = total_bmass$sal,
      y = total_bmass$um3
    ),
    size = 3, color = '#dbf6f6'
  )+
  geom_error_range(sal_pred, total_d, '#ff9191')+
  labs(x = "", y = "")+
  scale_y_log10() +
  annotation_logticks(sides = 'l')+
  theme_pubclean()+
  theme(
    axis.text = element_text(size = 20, color = 'white'),
    panel.grid = element_line(color = 'grey50'),
    plot.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.background = element_rect(fill = 'transparent')
  )
ggsave('./output/02b-tot_diat01-bmass_pred.pdf',
height = 6, width = 6, units = 'in')




# region \- sp changes  - - - - - - - - -
# size changes to solo pennate
sp_size = d_pred[[7]] |> 
  summarize_pred()

ggplot() +
  geom_point(
    data = indv_metrics |> 
      filter(group == 'solo_pennate'),
    aes(
      x = sal,
      y = um3
    )
  ) +
  geom_line(
    data = sp_size,
    aes(
      x = sal_pred,
      y = mean
    )
  ) +
  geom_ribbon(
    data = sp_size,
    aes(
      x = sal_pred,
      ymin = low,
      ymax = high
    ),
    alpha = 0.25
  ) +
  scale_y_log10()+
  annotation_logticks(sides = 'l') +
  theme_pubclean()


sp_count = n_pred[[7]] |> 
  summarize_pred()

ggplot() +
  geom_point(
    data = all_diat |> 
      filter(group == 'solo_pennate'),
    aes(
      x = sal,
      y = count
    )
  ) +
  geom_line(
    data = sp_count,
    aes(
      x = sal_pred,
      y = mean
    )
  ) +
  geom_ribbon(
    data = sp_count,
    aes(
      x = sal_pred,
      ymin = low,
      ymax = high
    ),
    alpha = 0.25
  ) +
  theme_pubclean()
