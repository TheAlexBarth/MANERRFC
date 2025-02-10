# temp plot sheet

all_conc <- conc |> 
  group_by(sample_id) |> 
  summarize(
    totalcount = sum(count),
    biomass = sum(um3_L),
    img_vol = sum(img_vol)
  )

ggplot() +
  geom_histogram(
    data = conc,
    aes(
      x = count,
      fill = taxa
    )
  ) +
  theme_bw()+
  theme(legend.position = 'none')

ggplot() +
  geom_density(
    data = all_conc,
    aes(x = biomass)
  ) +
  theme_bw()


ggplot() +
  geom_density(
    aes(
      x = all_conc$img_vol
    ),
    fill = 'black'
  )





  ggplot() +
    geom_point(
      aes(
        x = total_diat$sal,
        y = total_diat$biomass
      )
    )+
    geom_line(
      data = total_pred,
      aes(
        x = sal,
        y = mean
      )
    ) +
    geom_ribbon(
      data = total_pred,
      aes(
        x = sal,
        ymin = low,
        ymax = high
      ),
      alpha = 0.25
    ) +
    geom_line(
      data = total_lambda,
      aes(
        x = sal,
        y = mean
      ),
      color = 'red'
    ) +
    geom_ribbon(
      data = total_lambda,
      aes(
        x = sal,
        ymin = low,
        ymax = high
      ),
      fill = 'red',
      alpha = 0.25
    ) +
    labs(x = "", y = "")+
    theme_pubclean()    



# MARK: SALINITY 
ggplot(
  data = all_diat |> 
    group_by(year = year(date)) |> 
    summarize(
      s = mean(sal, na.rm = T),
      sig = sd(sal, na.rm = T)
    )
) +
  geom_point(
    aes(
      x = year,
      y = s
    ),
    size = 4
  ) +
  geom_line(
    aes(
      x = year,
      y = s
    )
  ) +
  geom_errorbar(
    aes(x = year, ymin = s - sig, ymax  = s+sig),
    width = 0.25
  ) +
  labs(x = 'Year', y = "Salinity") +
  theme_minimal()+
  theme(
    axis.title = element_text(face = 'bold', size = 16),
    axis.text = element_text(face = 'bold', size = 12)
  )



# MARK: TOTAL NPRED #############
  ggplot() +
    geom_point(
      aes(
        x = tobs_count$sal,
        y = tobs_count$count
      ),
      size = 3
    )+
    geom_line(
      data = total_npred,
      aes(
        x = sal,
        y = mean
      ),
      linewidth = 1.5
    ) +
    geom_ribbon(
      data = total_npred,
      aes(
        x = sal,
        ymin = low,
        ymax = high
      ),
      alpha = 0.25
    ) +
    # geom_line(
    #   data = total_lambda,
    #   aes(
    #     x = sal,
    #     y = mean
    #   ),
    #   color = 'red'
    # ) +
    # geom_ribbon(
    #   data = total_lambda,
    #   aes(
    #     x = sal,
    #     ymin = low,
    #     ymax = high
    #   ),
    #   fill = 'red',
    #   alpha = 0.25
    # ) +
    labs(x = "Salinity", y = "Diatom Count [#]")+
    scale_y_continuous(limits = c(0,18))+
    theme_pubclean() +
    theme(
      plot.background = element_blank(),
      axis.text = element_text(size = 12, face = 'bold'),
      axis.title = element_text(size = 16, face = 'bold')
    )



# MARK: BIOMASS TOTAL ###########################
ggplot() +
  geom_point(
    aes(
      x = total_bmass$sal,
      y = total_bmass$um3
    ),
    size = 3
  )+
  geom_line(
    data = total_d,
    aes(
      x = sal,
      y = mean
    ),
    linewidth = 1.5
  ) +
  geom_ribbon(
    data = total_d,
    aes(
      x = sal,
      ymin = low,
      ymax = high
    ),
    alpha = 0.25
  ) +
  scale_y_log10() +
  annotation_logticks(sides = 'l')+
  labs(
    x = "Salinity", 
    y = "Mean Individual Cell Biovolume [um3]"
  )+
  theme_pubclean() +
  theme(
    axis.title = element_text(face = 'bold', size = 16),
    axis.text = element_text(face = 'bold', size = 12)
  )


# MARK: TOTAL CONC ###############
ggplot() +
  geom_point(
    aes(
      x = total_diat$sal,
      y = total_diat$biomass
    ),
    size = 3
  )+
  geom_line(
    data = total_pred,
    aes(
      x = sal,
      y = mean
    ),
    linewidth = 1.5
  ) +
  geom_ribbon(
    data = total_pred,
    aes(
      x = sal,
      ymin = low,
      ymax = high
    ),
    alpha = 0.25
  ) +
  labs(
    x = "Salinity", 
    y = "Total Biovolume Concentration [um3/L]"
  )+
  geom_point(
    data = all_diat |> 
      filter(group == 'solo_pennate'),
    aes(
      x = sal,
      y = um3_L,
      color = gg_cbb_col(7)[7]
    ),
    size = 3
  ) +
  geom_line(
    aes(
      x = sal_pred,
      y = solo_pennate$mean
    ),
    linewidth = 1.5,
    color = gg_cbb_col(7)[7]
  ) +
  geom_ribbon(
    aes(
      x = sal_pred,
      ymin = solo_pennate$low,
      ymax = solo_pennate$high
    ),
    alpha = 0.25,
    fill = gg_cbb_col(7)[7]
  )+
  theme_pubclean() +
  scale_y_continuous(limits = c(0,4.15e10))+
  theme(
    axis.title = element_text(face = 'bold', size = 16),
    axis.text = element_text(face = 'bold', size = 12),
    legend.position = 'none'
  )




# MARK: PENNATE Conc #####################
ggplot() +
  geom_point(
    data = all_diat |> 
      filter(group == 'solo_pennate'),
    aes(
      x = sal,
      y = um3_L,
      color = gg_cbb_col(7)[7]
    ),
    size = 3
  ) +
  geom_line(
    aes(
      x = sal_pred,
      y = solo_pennate$mean
    ),
    linewidth = 1.5,
    color = gg_cbb_col(7)[7]
  ) +
  geom_ribbon(
    aes(
      x = sal_pred,
      ymin = solo_pennate$low,
      ymax = solo_pennate$high
    ),
    alpha = 0.05,
    fill = gg_cbb_col(7)[7]
  ) +
  labs(
    x = "Salinity", 
    y = "Total Biovolume Concentration [um3/L]"
  )+
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
    alpha = 0.5,
    fill = gg_cbb_col(8)[8]
  ) +
  theme_pubclean()+
  theme(legend.position = 'none') +
  theme(
    axis.title = element_text(face = 'bold', size = 16),
    axis.text = element_text(face = 'bold', size = 12)
  )
  



  ggplot() +
    geom_point(
      data = indv_metrics |> 
        filter(group == "solo_pennate"),
      aes(
        x = sal,
        y =um3
      ),
      size = 3,
      color = gg_cbb_col(7)[7]
    ) + 
    geom_line(
      data = pen_bmass,
      aes(
        x = sal,
        y = mean
      ),
      linewidth = 1.5,
      color = gg_cbb_col(7)[7]
    ) +
    geom_ribbon(
      data = pen_bmass,
      aes(
        x = sal,
        ymin = low,
        ymax = high
      ),
      fill = gg_cbb_col(7)[7],
      alpha = 0.25
    )+
    labs(x = "Salinity", y = "Cell Biovolume [um3]")+
    theme_pubclean()+
    theme(legend.position = 'none') +
    theme(
      axis.title = element_text(face = 'bold', size = 16),
      axis.text = element_text(face = 'bold', size = 12)
    )
      