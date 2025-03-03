rm(list = ls())
library(ggplot2)
library(dplyr)
library(rstan)
source('./R/utils.R')

# region \- load data ------------------
diat <- readRDS('./data/msi-diat_fit.rds')
wq = readRDS('./output/s02-temp-wq.rds')

# region \- water quality prep for ts ------------------

wq = wq |> 
    filter(sample_site %in% c('CE', 'CW')) |> 
    group_by(Date) |> 
    summarise(sal = mean(Sal), chl = mean(ChlFluor, na.rm = T))

wq = wq[!is.na(wq$sal),]


##########################################
# MARK: Post Pred fx---------------------
################################################

post_pred <- function(data, sal_range) {

  g_list <- list() #group list
  lambda <- list() # lambda list
  n_list <- list() # count list
  # loop through all levels of diatoms
  pb = txtProgressBar(0, length(levels(data$data$group)))
  setTxtProgressBar(pb, 0)
  for(l in 1:length(levels(data$data$group))) {
    #make lambda matrix (draws of posterior estimates for param)
    lambda[[l]] <- post_predict(
      rstan::extract(data$mod),
      list(
        sal = sal_range
      ),
      paste0(
        'exp(beta0[,',l,'] + beta1[,',l,'] * sal + log(0.1129333/1000))'
      ),
      n_draws = 2000
    )
    #from lambda mat - simulate counts
    n_list[[l]] <- matrix(
      rpois(length(lambda[[l]]), lambda = lambda[[l]]),
      nrow = nrow(lambda[[l]]),
      ncol = ncol(lambda[[l]])
    )
  }

  # between-group errors are v 
  err <- rstan::extract(data$mod, pars = 'v')$v
  eta <- list() # this is for biomass mean
  d_list <- list() # storage of d

  for(l in 1:length(levels(data$data$group))) {
    # draw group-specific mean (post param fitted estimates)
    eta[[l]] <- post_predict(
      rstan::extract(data$mod),
      list(
        sal = sal_range
      ),
      paste0(
        'a0[,',l,'] + a1[,',l,'] * sal'
      ),
      n_draws = 2000
    )
    # simulate from norm - biomass distribution
    d_list[[l]] <- apply(
      eta[[l]], 2, 
      function(x) rnorm(nrow(eta[[l]]), mean = x, sd = err)
    ) |> 
      exp()

    # predict biomass from d and n.
    g_list[[l]] = d_list[[l]] * (n_list[[l]]/mean(data$data$img_vol/1000))
    setTxtProgressBar(pb, l)
  }
  return(list(n = n_list, d = d_list, g = g_list))
}


########################################################
# MARK: Component Pred --------------------
################################################################

sal_pred <- seq(min(diat$data$sal), max(diat$data$sal), length.out = 1000)

# eh slower to have it but cleaner code
components <- post_pred(diat, sal_pred)

# region \- count pred -----------
total_n <- diat$data |> 
  group_by(sample_id, date) |> 
  summarize(
    total_n = sum(count),
    sal = mean(sal, na.rm = T)
  )

total_n_pred <- components$n |> sum_list(sp_range = sal_pred)

# region \-\- n count plot --------------

ggplot() +
  geom_point(
    data = total_n,
    aes(
      x = sal,
      y = total_n
    ),
    size = 3, color = '#f2dbf6'
  )+
  geom_error_range(sal_pred, total_n_pred, color = '#cc91ff') +
  labs(x = "", y = '')+
  theme_minimal()+
  theme(
    axis.text = element_text(size = 28, color = 'white'),
    panel.grid = element_line(color = 'grey50'),
    plot.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.background = element_rect(fill = 'transparent')
  )

ggsave('./output/msi/diatoms_n_comp.pdf',
height = 6, width = 6, units = 'in', dpi = 600
)

# region \- cell pred ---------------

total_d <- diat$indv |> 
  group_by(sample_id) |> 
  summarize(pgC = mean(pgC), sal = mean(sal, na.rm = T))

total_d_pred <- components$d |> 
  mean_list(sp_range = sal_pred)

total_d_pred_log <- total_d_pred[, names(total_d_pred)[-1]] |> 
  apply(c(1,2), log) |> 
  as.data.frame()

# region \-\- cell biomass plot -------

ggplot() +
  geom_point(
    data = total_d,
    aes(
      x = sal,
      y = log(pgC)
    ),
    size = 3, color = '#f2dbf6'
  )+
  geom_error_range(sal_pred, total_d_pred_log, color = '#cc91ff') +
  labs(x = "", y = '')+
  theme_minimal()+
  theme(
    axis.text = element_text(size = 28, color = 'white'),
    panel.grid = element_line(color = 'grey50'),
    plot.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.background = element_rect(fill = 'transparent')
  )

ggsave('./output/msi/diatoms_d_comp.pdf',
height = 6, width = 6, units = 'in', dpi = 600
)

# region \- full conc ---------------

total_g <- diat$data |> 
  group_by(sample_id) |> 
  summarize(pgC_L = sum(pgC_L), sal = mean(sal, na.rm = T))

total_g_pred <- components$g |> 
  sum_list(sal_pred)


# region \-\- conc plot --------------

ggplot() +
  geom_point(
    data = total_g,
    aes(
      x = sal,
      y = pgC_L
    ),
    size = 3, color = '#f2dbf6'
  )+
  geom_error_range(sal_pred, total_g_pred, color = '#cc91ff') +
  labs(x = "", y = '')+
  theme_minimal()+
  theme(
    axis.text = element_text(size = 28, color = 'white'),
    panel.grid = element_line(color = 'grey50'),
    plot.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.background = element_rect(fill = 'transparent')
  )

ggsave('./output/msi/diatoms_g_comp.pdf',
height = 6, width = 6, units = 'in', dpi = 600
)

  
####################################################
# MARK: TS Pred -------------------------
####################################################

g_ts_pred <- post_pred(diat, wq$sal)$g
    
total_diat <- diat$data |> 
  group_by(sample_id, date) |> 
  summarize(
    biomass = sum(pgC_L, na.rm = T),
    sal = mean(sal, na.rm = T)
  )


total_ts <- g_ts_pred |> 
  sum_list(sp_range = wq$sal)


total_plot_ts <- wq |> 
  left_join(
    total_ts,
    by = 'sal'
  )


# region \- total plot ---------------

ggplot() +
  geom_point(
    data = total_diat,
    aes(
      x = date,
      y = biomass
    ),
    size = 3, color = '#f2dbf6'
  )+
  geom_error_range(total_plot_ts$Date, total_plot_ts, '#cc91ff') +
  scale_x_date(
    limits = c(min(total_plot_ts$Date), max(total_plot_ts$Date))
  )+
  labs(x = "", y = '')+
  theme_minimal()+
  theme(
    axis.text = element_text(size = 28, color = 'white'),
    panel.grid = element_line(color = 'grey50'),
    plot.background = element_rect(fill = 'transparent', color = 'transparent'),
    panel.background = element_rect(fill = 'transparent')
  )

ggsave('./output/msi/diat-ts.pdf',
height = 6, width = 22, units = 'in')


# region \- diat density distribution--------

ggplot() +
  geom_histogram(
    data = diat$data,
    aes(
      x = log(pgC_L+1),
      group = group,
      fill = group
    ),
    alpha = 0.5
  )+
  labs(y = "Count", x = "log(pgC/L)", fill = "")+
  theme_minimal()+
    theme(
      axis.text = element_text(size = 28, color = 'black'),
      panel.grid = element_line(color = 'grey50'),
      plot.background = element_rect(fill = 'transparent', color = 'transparent'),
      panel.background = element_rect(fill = 'transparent')
    )

ggsave('./output/msi/diat-histo.png', width = 8, height = 4, units = 'in')
