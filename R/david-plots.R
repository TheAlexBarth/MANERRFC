### Script for generating David's plots ###

### Clear environment and load packages
rm(list = ls())
library(ggplot2)
library(dplyr)
library(lubridate)

### Read in merged data from module 02
full <- readRDS("./data/02-full_merged.rds")

####################################################

###Whole Community Plots

  #1. Lab Chlorophyll vs. Temperature
whole_CHLANvTemp <- ggplot(full$conc) +
  geom_point(
    aes(
      x = t,
      y = CHLA_N,
      color = sample_site
    )
  ) +
  geom_smooth(
    aes(
      x = t,
      y = CHLA_N,
      color = sample_site
    ),
    method = 'lm'
  )+
  theme_minimal()

whole_CHLANvTemp + ggtitle("Chlorophyll a Concentration vs. Water Temperature") +
  xlab("Average Monthly Water Temperature (°C)") + 
  ylab("Average Chlorophyll Fluorescence (RFU)")

ggsave("whole_CHLANvTemp.svg", plot = whole_CHLANvTemp, width = 7, height = 5)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

  #2. Lab Chlorophyll vs. salinity
whole_CHLANvSal <- ggplot(full$conc) +
  geom_point(
    aes(
      x = sal,
      y = CHLA_N,
      color = sample_site
    )
  ) +
  geom_smooth(
    aes(
      x = sal,
      y = CHLA_N,
      color = sample_site
    ),
    method = 'lm'
  )+
  theme_minimal()

whole_CHLANvSal + ggtitle("Chlorophyll a Concentration vs. Salinity") +
  xlab("Average Monthly Salinity (PSU)") + 
  ylab("Average Chlorophyll Fluorescence (RFU)")

ggsave("whole_CHLANvSal.svg", plot = whole_CHLANvSal, width = 7, height = 5)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

  #3. Lab Chlorophyll Biomass vs. Wind Speed
whole_CHLANvWindSpeed <- ggplot(full$conc) +
  geom_point(
    aes(
      x = windspeed,
      y = CHLA_N,
      color = sample_site
    )
  ) +
  geom_smooth(
    aes(
      x = windspeed,
      y = CHLA_N,
      color = sample_site
    ),
    method = 'lm'
  )+
  theme_minimal()

whole_CHLANvWindSpeed + ggtitle("Chlorophyll a Concentration vs. Wind Speed") +
  xlab("Average Monthly Wind Speed (m/s)") + 
  ylab("Average Chlorophyll Fluorescence (RFU)")

ggsave("whole_CHLANvWindSpeed.svg", plot = whole_CHLANvWindSpeed, width = 7, height = 5)
###########################################################################

### Diatom Biomass Plots

  #4. Diatom biomass vs Temperature 
diatom_Biomass_v_Temp <- ggplot(full$conc |> 
         filter(taxa %in% full$names$diatom) |> 
         group_by(id, sample_site) |> 
         summarize(
           pgC_L = sum(pgC_L, na.rm = T),
           t = unique(t)
         )
        ) +
          geom_point(
            aes(
              x = t,
              y = log(pgC_L+1),
              color = sample_site
            )
          ) +
          geom_smooth(
            aes(
              x = t,
              y = log(pgC_L+1),
              color = sample_site
            ),
            method = 'lm'
          )+ 
          theme_minimal()

diatom_Biomass_v_Temp + ggtitle("Diatom Biomass vs. Water Temperature") +
  xlab("Average Monthly Water Temperature (°C)") + 
  ylab("log Diatom Biomass Concentration (log(pgC/L)")

ggsave("diatom_Biomass_v_Temp.svg", plot = diatom_Biomass_v_Temp, width = 7, height = 5)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

  #5. Diatom biomass vs Salinity
diatom_Biomass_v_Sal <- ggplot(full$conc |> 
        filter(taxa %in% full$names$diatom) |> 
        group_by(id, sample_site) |> 
        summarize(
          pgC_L = sum(pgC_L, na.rm = T),
          sal = unique(sal)
        )
        ) +
          geom_point(
            aes(
              x = sal,
              y = log(pgC_L+1),
              color = sample_site
            )
          ) +
          geom_smooth(
            aes(
              x = sal,
              y = log(pgC_L+1),
              color = sample_site
            ),
            method = 'lm'
          )+ 
          theme_minimal()

diatom_Biomass_v_Sal + ggtitle("Diatom Biomass vs. Salinity") +
  xlab("Average Monthly Salinity (PSU)") + 
  ylab("log Diatom Biomass Concentration (log(pgC/L)")

ggsave("diatom_Biomass_v_Sal.svg", plot = diatom_Biomass_v_Sal, width = 7, height = 5)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

 #6. Diatom biomass vs Wind Speed 
diatom_Biomass_v_Wind <- ggplot(full$conc |> 
         filter(taxa %in% full$names$diatom) |> 
         group_by(id, sample_site) |> 
         summarize(
           pgC_L = sum(pgC_L, na.rm = T),
           windspeed = unique(windspeed)
         )
        ) +
          geom_point(
            aes(
              x = windspeed,
              y = log(pgC_L+1),
              color = sample_site
            )
          ) +
          geom_smooth(
            aes(
              x = windspeed,
              y = log(pgC_L+1),
              color = sample_site
            ),
            method = 'lm'
          )+ 
          theme_minimal()

diatom_Biomass_v_Wind + ggtitle("Diatom Biomass vs. Wind Speed") +
  xlab("Average Monthly Wind Speed (m/s)") + 
  ylab("log Diatom Biomass Concentration (log(pgC/L)")

ggsave("diatom_Biomass_v_Wind.svg", plot = diatom_Biomass_v_Wind, width = 7, height = 5)
###########################################################################

### Dinoflagellate Biomass Plots

#7. Dinoflagellate biomass vs Temperature 
dino_Biomass_v_Temp <- ggplot(full$conc |> 
          filter(taxa %in% full$names$dino) |> 
          group_by(id, sample_site) |> 
          summarize(
            pgC_L = sum(pgC_L, na.rm = T),
            t = unique(t)
          )
        ) +
          geom_point(
            aes(
              x = t,
              y = log(pgC_L+1),
              color = sample_site
            )
          ) +
          geom_smooth(
            aes(
              x = t,
              y = log(pgC_L+1),
              color = sample_site
            ),
            method = 'lm'
          )+ 
          theme_minimal()

dino_Biomass_v_Temp + ggtitle("Dinoflagellate Biomass vs. Water Temperature") +
  xlab("Average Monthly Water Temperature (°C)") + 
  ylab("log Dinoflagellate Biomass Concentration (log(pgC/L)")

ggsave("dino_Biomass_v_Temp.svg", plot = dino_Biomass_v_Temp, width = 7, height = 5)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#8. Dinoflagellate biomass vs Salinity
dino_Biomass_v_Sal <- ggplot(full$conc |> 
           filter(taxa %in% full$names$dino) |> 
           group_by(id, sample_site) |> 
           summarize(
             pgC_L = sum(pgC_L, na.rm = T),
             sal = unique(sal)
           )
        ) +
          geom_point(
            aes(
              x = sal,
              y = log(pgC_L+1),
              color = sample_site
            )
          ) +
          geom_smooth(
            aes(
              x = sal,
              y = log(pgC_L+1),
              color = sample_site
            ),
            method = 'lm'
          )+ 
          theme_minimal()

dino_Biomass_v_Sal + ggtitle("Dinoflagellate Biomass vs. Salinity") +
  xlab("Average Monthly Salinity (PSU)") + 
  ylab("log Dinoflagellate Biomass Concentration (log(pgC/L)")

ggsave("dino_Biomass_v_Sal.svg", plot = dino_Biomass_v_Sal, width = 7, height = 5)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#9. Dinoflagellate biomass vs Wind Speed 
dino_Biomass_v_Wind <- ggplot(full$conc |> 
          filter(taxa %in% full$names$dino) |> 
          group_by(id, sample_site) |> 
          summarize(
            pgC_L = sum(pgC_L, na.rm = T),
            windspeed = unique(windspeed)
          )
        ) +
          geom_point(
            aes(
              x = windspeed,
              y = log(pgC_L+1),
              color = sample_site
            )
          ) +
          geom_smooth(
            aes(
              x = windspeed,
              y = log(pgC_L+1),
              color = sample_site
            ),
            method = 'lm'
          )+ 
          theme_minimal()

dino_Biomass_v_Wind + ggtitle("Dinoflagellate Biomass vs. Wind Speed") +
  xlab("Average Monthly Wind Speed (m/s)") + 
  ylab("log Dinoflagellate Biomass Concentration (log(pgC/L)")

ggsave("dino_Biomass_v_Wind.svg", plot = dino_Biomass_v_Wind, width = 7, height = 5)
###########################################################################

### Ciliate Biomass Plots

#10. Ciliate biomass vs Temperature 
cili_Biomass_v_Temp <- ggplot(full$conc |> 
        filter(taxa %in% full$names$mz) |> 
        group_by(id, sample_site) |> 
        summarize(
          pgC_L = sum(pgC_L, na.rm = T),
          t = unique(t)
        )
      ) +
        geom_point(
          aes(
            x = t,
            y = log(pgC_L+1),
            color = sample_site
          )
        ) +
        geom_smooth(
          aes(
            x = t,
            y = log(pgC_L+1),
            color = sample_site
          ),
          method = 'lm'
        )+ 
        theme_minimal()

cili_Biomass_v_Temp + ggtitle("Ciliate Biomass vs. Water Temperature") +
  xlab("Average Monthly Water Temperature (°C)") + 
  ylab("log Ciliate Biomass Concentration (log(pgC/L)")

ggsave("cili_Biomass_v_Temp.svg", plot = cili_Biomass_v_Temp, width = 7, height = 5)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#11. Ciliate biomass vs Salinity
cili_Biomass_v_Sal <- ggplot(full$conc |> 
         filter(taxa %in% full$names$mz) |> 
         group_by(id, sample_site) |> 
         summarize(
           pgC_L = sum(pgC_L, na.rm = T),
           sal = unique(sal)
         )
        ) +
          geom_point(
            aes(
              x = sal,
              y = log(pgC_L+1),
              color = sample_site
            )
          ) +
          geom_smooth(
            aes(
              x = sal,
              y = log(pgC_L+1),
              color = sample_site
            ),
            method = 'lm'
          )+ 
          theme_minimal()

cili_Biomass_v_Sal + ggtitle("Ciliate Biomass vs. Salinity") +
  xlab("Average Monthly Salinity (PSU)") + 
  ylab("log Ciliate Biomass Concentration (log(pgC/L)")

ggsave("cili_Biomass_v_Sal.svg", plot = cili_Biomass_v_Sal, width = 7, height = 5)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#12. Ciliate biomass vs Wind Speed 
cili_Biomass_v_Wind <- ggplot(full$conc |> 
        filter(taxa %in% full$names$mz) |> 
        group_by(id, sample_site) |> 
        summarize(
          pgC_L = sum(pgC_L, na.rm = T),
          windspeed = unique(windspeed)
        )
      ) +
        geom_point(
          aes(
            x = windspeed,
            y = log(pgC_L+1),
            color = sample_site
          )
        ) +
        geom_smooth(
          aes(
            x = windspeed,
            y = log(pgC_L+1),
            color = sample_site
          ),
          method = 'lm'
        )+ 
        theme_minimal()

cili_Biomass_v_Wind + ggtitle("Ciliate Biomass vs. Wind Speed") +
  xlab("Average Monthly Wind Speed (m/s)") + 
  ylab("log Ciliate Biomass Concentration (log(pgC/L)")

ggsave("cili_Biomass_v_Wind.svg", plot = cili_Biomass_v_Wind, width = 7, height = 5)
###########################################################################
