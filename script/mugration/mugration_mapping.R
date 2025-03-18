# Install required libraries
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("maps", quietly = TRUE)) install.packages("maps")
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
if (!requireNamespace("tidygeocoder", quietly = TRUE)) install.packages("tidygeocoder")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

# Ensure required libraries are loaded
library(ggplot2)
library(maps)
library(ggrepel)
library(tidygeocoder)
library(dplyr)

# Load the migration data
data <- read.csv("annotated_tree_events.csv")
colnames(data) <- c("Index", "EventTime", "Origin", "Destination")

# Extract unique locations from Origin and Destination
locations <- unique(c(data$Origin, data$Destination))

# Geocode the unique locations using OpenStreetMap
locations_df <- data.frame(Location = locations) %>%
  geocode(Location, method = "osm") # Ensure internet access for this step

# Merge geocoded latitude and longitude for Origin
data <- merge(data, locations_df, by.x = "Origin", by.y = "Location", all.x = TRUE)
colnames(data)[which(colnames(data) %in% c("lat", "long"))] <- c("lat_Origin", "long_Origin")

# Merge geocoded latitude and longitude for Destination
data <- merge(data, locations_df, by.x = "Destination", by.y = "Location", all.x = TRUE)
colnames(data)[which(colnames(data) %in% c("lat", "long"))] <- c("lat_Destination", "long_Destination")

# Convert EventTime to a usable date format (assumes decimal years)
data$Date <- as.Date(as.POSIXct((data$EventTime - 1970) * 365.25 * 24 * 60 * 60, origin = "1970-01-01"))

# Sort by Date and filter the first 100 movements
data <- data %>%
  arrange(Date) %>%  # Sort by Date
  slice(-1) %>%      # Remove the first row (Unknown origin location)
  slice_head(n = 100) # Select the first 100 rows

# Sort date in descending order for visualisation (showing origings at the forefront)
data <- data %>%
  arrange(desc(Date))

# Load world map data
world_map <- map_data("world") # Contains 'group' and 'region' columns


## For visualisation

# Generate 5 evenly spaced dates based on the range of the Date field
dates <- seq(min(data$Date, na.rm = TRUE), max(data$Date, na.rm = TRUE), length.out = 5)

# Define custom color palette
color_palette <- scale_color_gradientn(
#  colours = c('red4', '#DB0201', 'darkorange', '#A0CBAD', '#8FB1BE', 'royalblue4'),
  colours = c(alpha('red4'), alpha('#DB0201', 0.9), alpha('darkorange', 0.8), alpha('#A0CBAD', 0.7), alpha('#8FB1BE', 0.6), alpha('royalblue4', 0.5)),
  breaks = dates,
  labels = format(dates, "%Y-%m"),
  name = 'Inferred date of dispersal'
)

# Plot migration events on the map
mugration_map <- ggplot() +
  theme_void() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "lightgrey", color = "white", linewidth = 0.2) +
  geom_curve(data = data,
             aes(x = long_Origin, y = lat_Origin, 
                 xend = long_Destination, yend = lat_Destination,
                 color = Date),
             curvature = 0.2, linewidth = 0.5, arrow = arrow(type = "closed", length = unit(0.1, "cm"))) +
  color_palette +
  #  labs(title = "Omicron BA.*") +
  #       subtitle = "First 100 movements", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(
    legend.position = 'top',
    legend.direction = 'horizontal',
    plot.title = element_text(family = "Helvetica"),
    legend.title = element_text(size = 12, hjust = 0.5),
    legend.text = element_text(size = 12),
    panel.grid = element_blank(), # Remove gridlines
    axis.title = element_blank(), # Remove axis titles
    axis.text = element_blank(),  # Remove axis text
    axis.ticks = element_blank()  # Remove axis ticks
  ) +
  guides(colour = guide_colourbar(barwidth = 25, barheight = 0.4, title.position = 'top', title.hjust = 1, ticks.colour = "white",
                                  ticks.linewidth = 0.5), size='none') +
  coord_fixed(ylim = c(-50, 90)) # Exclude Antarctica)

# Display the map
print(mugration_map)

# Save mugration map as a pdf
ggsave("ba*_mugration_map.pdf", plot = mugration_map, device = "pdf", width = 5.6, height = 3, dpi = 300)