rm(list = ls())

library(sf)
library(ggplot2)
library(ggmap)
library(ggspatial)
library(rnaturalearth)
source('./R/utils.R')

swmp_sites <- read.csv('./data/swmp_station_meta.csv')
wind <- read.csv('./data/xx-wind_sites.csv')


box_path <- '~/Library/CloudStorage/Box-Box/TGCRC Plankton Food Webs/Data/Maps/MinorBays'
# all_coast <- st_read(
#   paste0(box_path, '/CF inshore grids/CFgrid_AllCoast_poly_fixed.shp')
# )

minor_bays <- st_read(
  paste0(box_path, '/MinorBays.shp')
)

ma <- minor_bays[minor_bays$MAJORBAY == 5,]

swmp_sites <- st_as_sf(swmp_sites, coords = c('lon','lat')) |> 
  st_set_crs(4081) |> 
  st_transform(crs = st_crs(ma))

wind <- st_as_sf(wind, coords = c('lon','lat')) |> 
  st_set_crs(4081) |> 
  st_transform(crs = st_crs(ma))


swmp_sites$station_abbv <- toupper(swmp_sites$station_abbv)

ggplot() +
  annotation_map_tile(type = "osm", zoom = 10) +  # satellite imagery
  geom_sf(
    data = swmp_sites,
    aes(
      color = station_abbv
    ),
    size = 5
  )+
  geom_sf(
    data = wind,
    aes(
      shape = site,
    ),
    color = 'black',
    size = 3
  )+
  scale_color_manual(values = site_cols)+
  scale_shape_manual(values = c(7,9,10))+
  theme_minimal(base_size = 8) +
  labs(shape = "Wind Gauge", color = "Sampling Site")


ggsave('./output/fig01-map.pdf', width = 85, height = 85, units = 'mm', dpi = 600)


us_states <- ne_states(country = "United States of America", returnclass = "sf")
texas <- us_states[us_states$name == "Texas", ]

ma_bbox <- st_bbox(ma) |> st_as_sfc()

ggplot() +
  geom_sf(data = texas, fill = "white", color = "black") +
  geom_sf(data = ma_bbox, fill = NA, color = "red", linewidth = 1) +
  theme_void() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

ggsave('./output/fig01b-texas.pdf', width = 85, height = 85, units = 'mm', dpi = 600)
