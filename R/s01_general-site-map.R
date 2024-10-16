###
# Maps of MANERR #####
###


rm(list = ls())
library(ggplot2)
library(sf)
library(maps)
library(mapdata)

library(ggmap)
library(ggplot2)
library(grid)

source('./R/utils.R')

# Create a data frame with the GPS points
points_df <- data.frame(
  Location = c('AB', 'CE', 'CW'),
  Latitude = c(27.97980, 28.13230, 28.08410),
  Longitude = c(-97.02870, -97.03440, -97.20090)
)

points_df <- points_df[order(points_df$Location,
                            c('AB','CE','CW','MB','SC')),]

# Define the center and zoom level for your map
center <- c(lon = mean(points_df$Longitude), lat = mean(points_df$Latitude))
zoom <- 10

# Get the satellite image from Google Maps
# Requires registering an API key
satellite_map <- get_googlemap(center = center, zoom = zoom, maptype = 'satellite')

# Create the plot
sat_map <- ggmap(satellite_map) +
  geom_point(data = points_df, aes(x = Longitude, y = Latitude,
                                   color = Location), size = 3) +
  geom_text(data = points_df, aes(x = Longitude, y = Latitude, label = Location), 
            hjust = 0, vjust = 1.5, color = "white", size = 3) +
  scale_color_manual(values = gg_cbb_col(5)) +
  labs(x = "", y = "") +
  theme_minimal()

# # Convert the data frame to a simple feature object
points_sf <- st_as_sf(points_df, coords = c("Longitude", "Latitude"), crs = 4326)
# 
# # Get the map data for Texas
texas_map_data <- map_data("state", region = "texas")
# 
# Define the bounding box for the inset
bbox <- st_bbox(points_sf)
bbox_sf <- st_as_sfc(bbox)


# Inset map showing the bounding box within Texas
inset_map <- ggplot() +
  geom_polygon(data = texas_map_data, aes(x = long, y = lat, group = group), fill = "grey80", color = "white") +
  geom_sf(data = bbox_sf, color = "transparent", size = .001, linetype = "dashed") +
  geom_rect(aes(xmin = -97.4,xmax = -96.6,
                ymin = 27.75, ymax = 28.4),
            color = 'black', fill = 'black') +
  theme_void() +
  theme(plot.background = element_blank()) + 
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.box.background = element_rect(fill='transparent',
                                             linetype = 0),
        panel.border = element_rect(fill = 'transparent', color = 'transparent'))

# pdf('./figs/00_study-site-map.pdf')

ggsave('./output/03_study-site.pdf',
       plot = sat_map +
         inset(
           grob = ggplotGrob(inset_map),
           xmin = -96.9, xmax = -96.6, ymin = 27.6, ymax = 27.9
         ),
       height = 3,
       width = 3)

# dev.off()