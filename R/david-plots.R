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

whole_CHLANvTemp <- whole_CHLANvTemp + ggtitle("Chlorophyll a Concentration vs. Water Temperature") +
  xlab("Average Monthly Water Temperature (°C)") + 
  ylab("Average Chlorophyll Fluorescence (RFU)")

ggsave("whole_CHLANvTemp.png", plot = whole_CHLANvTemp, width = 7, height = 5)
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

whole_CHLANvSal <- whole_CHLANvSal + ggtitle("Chlorophyll a Concentration vs. Salinity") +
  xlab("Average Monthly Salinity (PSU)") + 
  ylab("Average Chlorophyll Fluorescence (RFU)")

ggsave("whole_CHLANvSal.png", plot = whole_CHLANvSal, width = 7, height = 5)
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

whole_CHLANvWindSpeed <- whole_CHLANvWindSpeed + ggtitle("Chlorophyll a Concentration vs. Wind Speed") +
  xlab("Average Monthly Wind Speed (m/s)") + 
  ylab("Average Chlorophyll Fluorescence (RFU)")

ggsave("whole_CHLANvWindSpeed.png", plot = whole_CHLANvWindSpeed, width = 7, height = 5)
ggsave("whole_CHLANvWindSpeed.svg", plot = whole_CHLANvWindSpeed, width = 7, height = 5)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#4. Lab Chlorophyll Biomass vs. Turbidity
whole_CHLAN_v_Turbid <- ggplot(full$conc) +
  geom_point(
    aes(
      x = Turb,
      y = CHLA_N,
      color = sample_site
    )
  ) +
  geom_smooth(
    aes(
      x = Turb,
      y = CHLA_N,
      color = sample_site
    ),
    method = 'lm'
  )+
  theme_minimal()

whole_CHLAN_v_Turbid <- whole_CHLAN_v_Turbid + ggtitle("Chlorophyll a Concentration vs. Turbidity") +
  xlab("Average Monthly Turbidity (FNU/NTU)") + 
  ylab("Average Chlorophyll Fluorescence (RFU)")

ggsave("whole_CHLAN_v_Turbid.png", plot = whole_CHLAN_v_Turbid, width = 7, height = 5)
ggsave("whole_CHLAN_v_Turbid.svg", plot = whole_CHLAN_v_Turbid, width = 7, height = 5)
###########################################################################

### Diatom Biomass Plots

  #5. Diatom biomass vs Temperature 
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

diatom_Biomass_v_Temp <- diatom_Biomass_v_Temp + ggtitle("Diatom Biomass Concentration vs. Water Temperature") +
  xlab("Average Monthly Water Temperature (°C)") + 
  ylab("log Diatom Biomass Concentration (log(pgC/L)")

ggsave("diatom_Biomass_v_Temp.png", plot = diatom_Biomass_v_Temp, width = 7, height = 5)
ggsave("diatom_Biomass_v_Temp.svg", plot = diatom_Biomass_v_Temp, width = 7, height = 5)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

  #6. Diatom biomass vs Salinity
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

diatom_Biomass_v_Sal <- diatom_Biomass_v_Sal + ggtitle("Diatom Biomass Concentration vs. Salinity") +
  xlab("Average Monthly Salinity (PSU)") + 
  ylab("log Diatom Biomass Concentration (log(pgC/L)")

ggsave("diatom_Biomass_v_Sal.png", plot = diatom_Biomass_v_Sal, width = 7, height = 5)
ggsave("diatom_Biomass_v_Sal.svg", plot = diatom_Biomass_v_Sal, width = 7, height = 5)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

 #7. Diatom biomass vs Wind Speed 
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
            )
          )+ 
          theme_minimal()

diatom_Biomass_v_Wind <- diatom_Biomass_v_Wind + ggtitle("Diatom Biomass Concentration vs. Wind Speed") +
  xlab("Average Monthly Wind Speed (m/s)") + 
  ylab("log Diatom Biomass Concentration (log(pgC/L)")

diatom_Biomass_v_Wind

ggsave("diatom_Biomass_v_Wind.png", plot = diatom_Biomass_v_Wind, width = 7, height = 5)
ggsave("diatom_Biomass_v_Wind.svg", plot = diatom_Biomass_v_Wind, width = 7, height = 5)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#8. Diatom Biomass vs. Turbidity
diatom_Biomass_v_Turbid <- ggplot(full$conc |> 
          filter(taxa %in% full$names$diatom) |> 
          group_by(id, sample_site) |> 
          summarize(
            pgC_L = sum(pgC_L, na.rm = T),
            Turb = unique(Turb)
          )
        ) +
          geom_point(
            aes(
              x = Turb,
              y = log(pgC_L+1),
              color = sample_site
            )
          ) +
          geom_smooth(
            aes(
              x = Turb,
              y = log(pgC_L+1),
              color = sample_site
            ),
            method = 'lm'
          )+ 
          theme_minimal()
#diatom_Biomass_v_Turbid

diatom_Biomass_v_Turbid <- diatom_Biomass_v_Turbid + ggtitle("Diatom Biomass Concentration vs. Turbidity") +
  xlab("Average Monthly Turbidity (FNU/NTU)") + 
  ylab("log Diatom Biomass Concentration (log(pgC/L)")

ggsave("diatom_Biomass_v_Turbid.png", plot = diatom_Biomass_v_Turbid, width = 7, height = 5)
ggsave("diatom_Biomass_v_Turbid.svg", plot = diatom_Biomass_v_Turbid, width = 7, height = 5)
###########################################################################

### Dinoflagellate Biomass Plots

#9. Dinoflagellate biomass vs Temperature 
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

dino_Biomass_v_Temp <- dino_Biomass_v_Temp + ggtitle("Dinoflagellate Biomass Concentration vs. Water Temperature") +
  xlab("Average Monthly Water Temperature (°C)") + 
  ylab("log Dinoflagellate Biomass Concentration (log(pgC/L)")

ggsave("dino_Biomass_v_Temp.png", plot = dino_Biomass_v_Temp, width = 7, height = 5)
ggsave("dino_Biomass_v_Temp.svg", plot = dino_Biomass_v_Temp, width = 7, height = 5)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#10. Dinoflagellate biomass vs Salinity
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

dino_Biomass_v_Sal <- dino_Biomass_v_Sal + ggtitle("Dinoflagellate Biomass Concentration vs. Salinity") +
  xlab("Average Monthly Salinity (PSU)") + 
  ylab("log Dinoflagellate Biomass Concentration (log(pgC/L)")

ggsave("dino_Biomass_v_Sal.png", plot = dino_Biomass_v_Sal, width = 7, height = 5)
ggsave("dino_Biomass_v_Sal.svg", plot = dino_Biomass_v_Sal, width = 7, height = 5)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#11. Dinoflagellate biomass vs Wind Speed 
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

dino_Biomass_v_Wind <- dino_Biomass_v_Wind + ggtitle("Dinoflagellate Biomass Concentration vs. Wind Speed") +
  xlab("Average Monthly Wind Speed (m/s)") + 
  ylab("log Dinoflagellate Biomass Concentration (log(pgC/L)")

ggsave("dino_Biomass_v_Wind.png", plot = dino_Biomass_v_Wind, width = 7, height = 5)
ggsave("dino_Biomass_v_Wind.svg", plot = dino_Biomass_v_Wind, width = 7, height = 5)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#12. Dinoflagellate Biomass vs. Turbidity
dino_Biomass_v_Turbid <- ggplot(full$conc |> 
        filter(taxa %in% full$names$dino) |> 
        group_by(id, sample_site) |> 
        summarize(
          pgC_L = sum(pgC_L, na.rm = T),
          Turb = unique(Turb)
        )
      ) +
        geom_point(
          aes(
            x = Turb,
            y = log(pgC_L+1),
            color = sample_site
          )
        ) +
        geom_smooth(
          aes(
            x = Turb,
            y = log(pgC_L+1),
            color = sample_site
          ),
          method = 'lm'
        )+ 
        theme_minimal()
#dino_Biomass_v_Turbid

dino_Biomass_v_Turbid <- dino_Biomass_v_Turbid + ggtitle("Dinoflagellate Biomass Concentration vs. Turbidity") +
  xlab("Average Monthly Turbidity (FNU/NTU)") + 
  ylab("log Dinoflagellate Biomass Concentration (log(pgC/L)")

ggsave("dino_Biomass_v_Turbid.png", plot = dino_Biomass_v_Turbid, width = 7, height = 5)
ggsave("dino_Biomass_v_Turbid.svg", plot = dino_Biomass_v_Turbid, width = 7, height = 5)
###########################################################################

### Ciliate Biomass Plots

#13. Ciliate biomass vs Temperature 
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

cili_Biomass_v_Temp <- cili_Biomass_v_Temp + ggtitle("Ciliate Biomass Concentration vs. Water Temperature") +
  xlab("Average Monthly Water Temperature (°C)") + 
  ylab("log Ciliate Biomass Concentration (log(pgC/L)")

ggsave("cili_Biomass_v_Temp.png", plot = cili_Biomass_v_Temp, width = 7, height = 5)
ggsave("cili_Biomass_v_Temp.svg", plot = cili_Biomass_v_Temp, width = 7, height = 5)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#14. Ciliate biomass vs Salinity
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

cili_Biomass_v_Sal <- cili_Biomass_v_Sal + ggtitle("Ciliate Biomass Concentration vs. Salinity") +
  xlab("Average Monthly Salinity (PSU)") + 
  ylab("log Ciliate Biomass Concentration (log(pgC/L)")

ggsave("cili_Biomass_v_Sal.png", plot = cili_Biomass_v_Sal, width = 7, height = 5)
ggsave("cili_Biomass_v_Sal.svg", plot = cili_Biomass_v_Sal, width = 7, height = 5)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#15. Ciliate biomass vs Wind Speed 
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

cili_Biomass_v_Wind <- cili_Biomass_v_Wind + ggtitle("Ciliate Biomass Concentration vs. Wind Speed") +
  xlab("Average Monthly Wind Speed (m/s)") + 
  ylab("log Ciliate Biomass Concentration (log(pgC/L)")

ggsave("cili_Biomass_v_Wind.png", plot = cili_Biomass_v_Wind, width = 7, height = 5)
ggsave("cili_Biomass_v_Wind.svg", plot = cili_Biomass_v_Wind, width = 7, height = 5)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#16. Ciliate Biomass vs. Turbidity
cili_Biomass_v_Turbid <- ggplot(full$conc |> 
          filter(taxa %in% full$names$mz) |> 
          group_by(id, sample_site) |> 
          summarize(
            pgC_L = sum(pgC_L, na.rm = T),
            Turb = unique(Turb)
          )
        ) +
          geom_point(
            aes(
              x = Turb,
              y = log(pgC_L+1),
              color = sample_site
            )
          ) +
          geom_smooth(
            aes(
              x = Turb,
              y = log(pgC_L+1),
              color = sample_site
            ),
            method = 'lm'
          )+ 
          theme_minimal()


cili_Biomass_v_Turbid <- cili_Biomass_v_Turbid + ggtitle("Ciliate Biomass Concentration vs. Turbidity") +
  xlab("Average Monthly Turbidity (FNU/NTU)") + 
  ylab("log Ciliate Biomass Concentration (log(pgC/L)")

ggsave("cili_Biomass_v_Turbid.png", plot = cili_Biomass_v_Turbid, width = 7, height = 5)
ggsave("cili_Biomass_v_Turbid.svg", plot = cili_Biomass_v_Turbid, width = 7, height = 5)
###########################################################################

#17 Stacked Bar Plot for whole estuary, taxa colored by diatom/dino/mz, grouped to year

#create dataset

full$conc <- full$conc |> 
  mutate(
    category = case_when(
      taxa %in% full$names$diatom ~ "Diatoms",
      taxa %in% full$names$dino ~ "Dinos",
      taxa %in% full$names$mz ~ "MZ",
      TRUE ~ NA_character_  # Exclude taxa that don’t match
    )
  ) |> 
  filter(!is.na(category))  # Remove non-matching taxa
# Summarize biomass by category and time
biomass_summary <- full$conc |> 
  group_by(yearmo, category) |> 
  summarize(
    total_pgC_L = sum(pgC_L, na.rm = TRUE),
    .groups = "drop"
  )
# Create stacked bar plot
ggplot(biomass_summary, aes(x = yearmo, y = total_pgC_L, fill = category)) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked bars
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +  # Format x-axis
  labs(title = "Stacked Biomass Time Series",
       x = "Time (Year-Month)", 
       y = "Biomass (pgC/L)", 
       fill = "Taxa Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#18 Stacked Bar Plot for Aransas Bay, taxa colored by diatom/dino/mz, grouped to year

#create dataset

AB_Stack_Filter <- full$conc |> 
  filter(sample_site == "AB") |>  # Filter for sample site 'AB'
  mutate(
    category = case_when(
      taxa %in% full$names$diatom ~ "Diatoms",
      taxa %in% full$names$dino ~ "Dinos",
      taxa %in% full$names$mz ~ "MZ",
      TRUE ~ NA_character_  # Exclude taxa that don’t match
    )
  ) |> 
  filter(!is.na(category))  # Remove non-matching taxa
# Summarize biomass by category and time
biomass_summary <- AB_Stack_Filter |> 
  group_by(yearmo, category) |> 
  summarize(
    total_pgC_L = sum(pgC_L, na.rm = TRUE),
    .groups = "drop"
  )
# Create stacked bar plot
ggplot(biomass_summary, aes(x = yearmo, y = total_pgC_L, fill = category)) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked bars
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +  # Format x-axis
  labs(title = "Stacked Biomass Time Series (Sample Site: AB)",
       x = "Time (Year-Month)", 
       y = "Biomass (pgC/L)", 
       fill = "Taxa Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#19 Stacked Bar Plot for Copano West, taxa colored by diatom/dino/mz, grouped to year

#create dataset

CW_Stack_Filter <- full$conc |> 
  filter(sample_site == "CW") |>  # Filter for sample site 'CW'
  mutate(
    category = case_when(
      taxa %in% full$names$diatom ~ "Diatoms",
      taxa %in% full$names$dino ~ "Dinos",
      taxa %in% full$names$mz ~ "MZ",
      TRUE ~ NA_character_  # Exclude taxa that don’t match
    )
  ) |> 
  filter(!is.na(category))  # Remove non-matching taxa
# Summarize biomass by category and time
biomass_summary <- CW_Stack_Filter |> 
  group_by(yearmo, category) |> 
  summarize(
    total_pgC_L = sum(pgC_L, na.rm = TRUE),
    .groups = "drop"
  )
# Create stacked bar plot
CW_Stack <- ggplot(biomass_summary, aes(x = yearmo, y = total_pgC_L, fill = category)) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked bars
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +  # Format x-axis
  labs(title = "Stacked Biomass Time Series (Sample Site: CW)",
       x = "Time (Year-Month)", 
       y = "% Biomass", 
       fill = "Taxa Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

ggsave("CW_Stack.png", plot = CW_Stack, width = 13, height = 2.25)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#20 Stacked Bar Plot for Copano East, taxa colored by diatom/dino/mz, grouped to year

#create dataset

CE_Stack_Filter <- full$conc |> 
  filter(sample_site == "CE") |>  # Filter for sample site 'CE'
  mutate(
    category = case_when(
      taxa %in% full$names$diatom ~ "Diatoms",
      taxa %in% full$names$dino ~ "Dinos",
      taxa %in% full$names$mz ~ "MZ",
      TRUE ~ NA_character_  # Exclude taxa that don’t match
    )
  ) |> 
  filter(!is.na(category))  # Remove non-matching taxa
# Summarize biomass by category and time
biomass_summary <- CE_Stack_Filter |> 
  group_by(yearmo, category) |> 
  summarize(
    total_pgC_L = sum(pgC_L, na.rm = TRUE),
    .groups = "drop"
  )
# Create stacked bar plot
ggplot(biomass_summary, aes(x = yearmo, y = total_pgC_L, fill = category)) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked bars
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +  # Format x-axis
  labs(title = "Stacked Biomass Time Series (Sample Site: CE)",
       x = "Time (Year-Month)", 
       y = "Biomass (pgC/L)", 
       fill = "Taxa Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#21 Stacked Bar Plot for Mesquite Bay, taxa colored by diatom/dino/mz, grouped to year

#create dataset

MB_Stack_Filter <- full$conc |> 
  filter(sample_site == "MB") |>  # Filter for sample site 'MB'
  mutate(
    category = case_when(
      taxa %in% full$names$diatom ~ "Diatoms",
      taxa %in% full$names$dino ~ "Dinos",
      taxa %in% full$names$mz ~ "MZ",
      TRUE ~ NA_character_  # Exclude taxa that don’t match
    )
  ) |> 
  filter(!is.na(category))  # Remove non-matching taxa
# Summarize biomass by category and time
biomass_summary <- MB_Stack_Filter |> 
  group_by(yearmo, category) |> 
  summarize(
    total_pgC_L = sum(pgC_L, na.rm = TRUE),
    .groups = "drop"
  )
# Create stacked bar plot
ggplot(biomass_summary, aes(x = yearmo, y = total_pgC_L, fill = category)) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked bars
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +  # Format x-axis
  labs(title = "Stacked Biomass Time Series (Sample Site: MB)",
       x = "Time (Year-Month)", 
       y = "Biomass (pgC/L)", 
       fill = "Taxa Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#22 Stacked Bar Plot for Ship Channel, taxa colored by diatom/dino/mz, grouped to year

#create dataset

SC_Stack_Filter <- full$conc |> 
  filter(sample_site == "SC") |>  # Filter for sample site 'SC'
  mutate(
    category = case_when(
      taxa %in% full$names$diatom ~ "Diatoms",
      taxa %in% full$names$dino ~ "Dinos",
      taxa %in% full$names$mz ~ "MZ",
      TRUE ~ NA_character_  # Exclude taxa that don’t match
    )
  ) |> 
  filter(!is.na(category))  # Remove non-matching taxa
# Summarize biomass by category and time
biomass_summary <- SC_Stack_Filter |> 
  group_by(month, category) |> 
  summarize(
    total_pgC_L = mean(pgC_L, na.rm = TRUE),
    .groups = "drop"
  )
# Create stacked bar plot
ggplot(biomass_summary, aes(x = month, y = total_pgC_L, fill = category)) +
  geom_bar(stat = "identity", position = "dodge") +  # Stacked bars
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +  # Format x-axis
  labs(title = "Stacked Biomass Time Series (Sample Site: SC)",
       x = "Time (Year-Month)", 
       y = "Biomass (pgC/L)", 
       fill = "Taxa Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#23 Stacked Bar Plot for Copano West, taxa colored by diatom/dino/mz, grouped by month, relative biomass abundance

#create dataset

CW_Stack_Filter <- full$conc |> 
  filter(sample_site == "CW") |>  # Filter for sample site 'CW'
  mutate(
    category = case_when(
      taxa %in% full$names$diatom ~ "Diatoms",
      taxa %in% full$names$dino ~ "Dinos",
      taxa %in% full$names$mz ~ "MZ",
      TRUE ~ NA_character_  # Exclude taxa that don’t match
    )
  ) |> 
  filter(!is.na(category))  # Remove non-matching taxa
# Summarize biomass by category and time
biomass_summary <- CW_Stack_Filter |> 
  group_by(yearmo, category) |> 
  summarize(
    total_pgC_L = sum(pgC_L, na.rm = TRUE),
    .groups = "drop"
  )
# Create stacked bar plot
CW_Stack <- ggplot(biomass_summary, aes(x = yearmo, y = total_pgC_L, fill = category)) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked bars
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +  # Format x-axis
  labs(title = "River Endmember Plankton Relative Abundance (Site CW)",
       x = "Time (Year-Month)", 
       y = "% Biomass", 
       fill = "Taxa Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

ggsave("CW_Stack.png", plot = CW_Stack, width = 13, height = 2.25)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#24 Stacked Bar Plot for Mesquite Bay, taxa colored by diatom/dino/mz, grouped by month, relative biomass abundance

#create dataset

MB_Stack_Filter <- full$conc |> 
  filter(sample_site == "MB") |>  # Filter for sample site 'MB'
  mutate(
    category = case_when(
      taxa %in% full$names$diatom ~ "Diatoms",
      taxa %in% full$names$dino ~ "Dinos",
      taxa %in% full$names$mz ~ "MZ",
      TRUE ~ NA_character_  # Exclude taxa that don’t match
    )
  ) |> 
  filter(!is.na(category))  # Remove non-matching taxa
# Summarize biomass by category and time
biomass_summary <- MB_Stack_Filter |> 
  group_by(yearmo, category) |> 
  summarize(
    total_pgC_L = sum(pgC_L, na.rm = TRUE),
    .groups = "drop"
  )
# Create stacked bar plot
MB_Stack <- ggplot(biomass_summary, aes(x = yearmo, y = total_pgC_L, fill = category)) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked bars
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +  # Format x-axis
  labs(title = "Intermediate Salinity Plankton Relative Abundance (Site MB)",
       x = "Time (Year-Month)", 
       y = "% Biomass", 
       fill = "Taxa Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

ggsave("MB_Stack.png", plot = MB_Stack, width = 13, height = 2.25)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#25 Stacked Bar Plot for Ship Channel, taxa colored by diatom/dino/mz, grouped by month, relative biomass abundance

#create dataset

SC_Stack_Filter <- full$conc |> 
  filter(sample_site == "SC") |>  # Filter for sample site 'SC'
  mutate(
    category = case_when(
      taxa %in% full$names$diatom ~ "Diatoms",
      taxa %in% full$names$dino ~ "Dinos",
      taxa %in% full$names$mz ~ "MZ",
      TRUE ~ NA_character_  # Exclude taxa that don’t match
    )
  ) |> 
  filter(!is.na(category))  # Remove non-matching taxa
# Summarize biomass by category and time
biomass_summary <- SC_Stack_Filter |> 
  group_by(yearmo, category) |> 
  summarize(
    total_pgC_L = sum(pgC_L, na.rm = TRUE),
    .groups = "drop"
  )
# Create stacked bar plot
SC_Stack <- ggplot(biomass_summary, aes(x = yearmo, y = total_pgC_L, fill = category)) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked bars
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +  # Format x-axis
  labs(title = "Marine Endmember Plankton Relative Abundance (Site SC)",
       x = "Time (Year-Month)", 
       y = "% Biomass", 
       fill = "Taxa Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

ggsave("SC_Stack.png", plot = SC_Stack, width = 13, height = 2.25)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#26 Stacked Bar Plot for Copano West, taxa colored by diatom/dino/mz, grouped by month, relative biomass abundance

#create dataset

CW_Stack_Filter <- full$conc |> 
  filter(sample_site == "CW") |>  # Filter for sample site 'CW'
  mutate(
    category = case_when(
      taxa %in% full$names$diatom ~ "Diatoms",
      taxa %in% full$names$dino ~ "Dinos",
      taxa %in% full$names$mz ~ "MZ",
      TRUE ~ NA_character_  # Exclude taxa that don’t match
    )
  ) |> 
  filter(!is.na(category))  # Remove non-matching taxa
# Summarize biomass by category and time
biomass_summary <- CW_Stack_Filter |> 
  group_by(month, category) |> 
  summarize(
    total_pgC_L = sum(pgC_L, na.rm = TRUE),
    .groups = "drop"
  )
# Create stacked bar plot
CW_Dodge_Month <- ggplot(biomass_summary, aes(x = month, y = total_pgC_L, fill = category)) +
  geom_bar(stat = "identity", position = "dodge") +  # Stacked bars
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +  # Format x-axis
  scale_y_continuous(limits = c(0, 1.0e+10), labels = scales::scientific) + 
  labs(title = "River Endmember Average Plankton Biomass (Site CW)",
       x = "Time (Month)", 
       y = "Biomass (PgC/L)", 
       fill = "Taxa Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

CW_Dodge_Month

ggsave("CW_Dodge_Month.png", plot = CW_Dodge_Month, width = 6, height = 6)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#27 Stacked Bar Plot for Mesquite Bay, taxa colored by diatom/dino/mz, grouped by month, relative biomass abundance

#create dataset

MB_Stack_Filter <- full$conc |> 
  filter(sample_site == "MB") |>  # Filter for sample site 'MB'
  mutate(
    category = case_when(
      taxa %in% full$names$diatom ~ "Diatoms",
      taxa %in% full$names$dino ~ "Dinos",
      taxa %in% full$names$mz ~ "MZ",
      TRUE ~ NA_character_  # Exclude taxa that don’t match
    )
  ) |> 
  filter(!is.na(category))  # Remove non-matching taxa
# Summarize biomass by category and time
biomass_summary <- MB_Stack_Filter |> 
  group_by(month, category) |> 
  summarize(
    total_pgC_L = sum(pgC_L, na.rm = TRUE),
    .groups = "drop"
  )
# Create stacked bar plot
MB_Dodge_Month <- ggplot(biomass_summary, aes(x = month, y = total_pgC_L, fill = category)) +
  geom_bar(stat = "identity", position = "dodge") +  # Stacked bars
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +  # Format x-axis
  scale_y_continuous(limits = c(0, 1.0e+10), labels = scales::scientific) + 
  labs(title = "Intermediate Salinity Average Plankton Biomass (Site MB)",
       x = "Time (Month)", 
       y = "Biomass (PgC/L)", 
       fill = "Taxa Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

MB_Dodge_Month

ggsave("MB_Dodge_Month.png", plot = MB_Dodge_Month, width = 6, height = 6)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#28 Stacked Bar Plot for Ship Channel, taxa colored by diatom/dino/mz, grouped by month, relative biomass abundance

#create dataset

SC_Stack_Filter <- full$conc |> 
  filter(sample_site == "SC") |>  # Filter for sample site 'SC'
  mutate(
    category = case_when(
      taxa %in% full$names$diatom ~ "Diatoms",
      taxa %in% full$names$dino ~ "Dinos",
      taxa %in% full$names$mz ~ "MZ",
      TRUE ~ NA_character_  # Exclude taxa that don’t match
    )
  ) |> 
  filter(!is.na(category))  # Remove non-matching taxa
# Summarize biomass by category and time
biomass_summary <- SC_Stack_Filter |> 
  group_by(month, category) |> 
  summarize(
    total_pgC_L = sum(pgC_L, na.rm = TRUE),
    .groups = "drop"
  )
# Create stacked bar plot
SC_Dodge_Month <- ggplot(biomass_summary, aes(x = month, y = total_pgC_L, fill = category)) +
  geom_bar(stat = "identity", position = "dodge") +  # Stacked bars
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +  # Format x-axis
  scale_y_continuous(limits = c(0, 1.0e+10), labels = scales::scientific) + 
  labs(title = "Marine Endmember Average Plankton Biomass (Site SC)",
       x = "Time (Month)", 
       y = "Biomass (PgC/L)", 
       fill = "Taxa Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

SC_Dodge_Month

ggsave("SC_Dodge_Month.png", plot = SC_Dodge_Month, width = 6, height = 6)
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_