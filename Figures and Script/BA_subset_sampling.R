### Omicron BA subsetting and sampling ###


## Libraries
library(dplyr)
library(data.table)
library("readxl")
library(lubridate)
library(osmdata)
library(ggplot2)
library(ggmap)
library(leaflet)


## Set working directory
setwd("")

## Read in metadata
metadata <- read_excel("")
nextclade <- read_excel("")


## Join tables
metadata_nextclade <- left_join(metadata, nextclade, 
                                by = c("gisaid_epi_isl" = "gisaid_epi_isl_nextclade"))

non_matched_nextclade <- anti_join(nextclade, metadata, by = c("gisaid_epi_isl_nextclade" = "gisaid_epi_isl"))

## Sort data by date
metadata_nextclade <- metadata_nextclade[order(metadata_nextclade$date), ]

## Filter BA sequences
BA1_like <- subset(metadata_nextclade, metadata_nextclade$pango_lineage %like% "BA.1")
BA1 <- subset(metadata_nextclade, grepl("^BA\\.1(\\.[0-9]+)*$", metadata_nextclade$pango_lineage))
#BA2 <- subset(metadata_nextclade, metadata_nextclade$pango_lineage=="BA.2")
BA2 <- subset(metadata_nextclade, grepl("^BA\\.2(\\.[0-9]+)*$", metadata_nextclade$pango_lineage))
#BA3 <- subset(metadata_nextclade, metadata_nextclade$pango_lineage=="BA.3")
BA3 <- subset(metadata_nextclade, grepl("^BA\\.3(\\.[0-9]+)*$", metadata_nextclade$pango_lineage))
#BA4 <- subset(metadata_nextclade, metadata_nextclade$pango_lineage=="BA.4")
BA4 <- subset(metadata_nextclade, grepl("^BA\\.4(\\.[0-9]+)*$", metadata_nextclade$pango_lineage))
#BA5 <- subset(metadata_nextclade, metadata_nextclade$pango_lineage=="BA.5")
BA5 <- subset(metadata_nextclade, grepl("^BA\\.5(\\.[0-9]+)*$", metadata_nextclade$pango_lineage))
#BA6 <- subset(metadata_nextclade, metadata_nextclade$pango_lineage=="BA.2.86")
BA6 <- subset(metadata_nextclade, grepl("^BA\\.2\\.86(\\.[0-9]+)*$", metadata_nextclade$pango_lineage))


## Filter first four months or until there are over 100 observations (dates are based off first occurrence in dataset)
# BA1
start_date_ba1 <- date("2021-10-12")
end_date_ba1 <- start_date_ba1 %m+% months(4)
BA1_4month <- BA1 %>%
  filter(date >= start_date_ba1 & date < end_date_ba1)
#BA2
start_date_ba2 <- date("2021-11-17")
end_date_ba2 <- start_date_ba2 %m+% months(4)
BA2_4month <- BA2 %>%
  filter(date >= start_date_ba2 & date < end_date_ba2)
#BA3
start_date_ba3 <- date("2021-11-18")
end_date_ba3 <- start_date_ba3 %m+% months(4)
BA3_4month <- BA3 %>%
  filter(date >= start_date_ba3 & date < end_date_ba3)
#BA4
start_date_ba4 <- date("2022-01-10")
end_date_ba4 <- start_date_ba4 %m+% months(4)
BA4_4month <- BA4 %>%
  filter(date >= start_date_ba4 & date < end_date_ba4)
#BA5
start_date_ba5 <- date("2022-01-11")
end_date_ba5 <- start_date_ba5 %m+% months(4)
BA5_4month <- BA5 %>%
  filter(date >= start_date_ba5 & date < end_date_ba5)
#BA6
start_date_ba6 <- date("2023-07-24")
end_date_ba6 <- start_date_ba6 %m+% months(4)
BA6_4month <- BA6 %>%
  filter(date >= start_date_ba6 & date < end_date_ba6)



## Random sample of respective BA sequences

# BA1
set.seed(101) # Set seed for reproducibility
BA1_first20 <- BA1_4month[1:20, ] # Select the first 20 rows after sorting
BA1_random80 <- BA1_4month[21:nrow(BA1_4month), ][sample(nrow(BA1_4month) - 20, 80, replace = FALSE), ] # Randomly sample 80 rows from the remainder
BA1_sample <- rbind(BA1_first20, BA1_random80) # Combine the two sets
BA1_sample <- BA1_sample[order(BA1_sample$date), ] # Sort by date

# BA2
set.seed(201)
BA2_first20 <- BA2_4month[1:20, ]
BA2_random80 <- BA2_4month[21:nrow(BA2_4month), ][sample(nrow(BA2_4month) - 20, 80, replace = FALSE), ]
BA2_sample <- rbind(BA2_first20, BA2_random80)
BA2_sample <- BA2_sample[order(BA2_sample$date), ]

# BA3
set.seed(301)
BA3_first20 <- BA3_4month[1:20, ]
BA3_random80 <- BA3_4month[21:nrow(BA3_4month), ][sample(nrow(BA3_4month) - 20, 80, replace = FALSE), ]
BA3_sample <- rbind(BA3_first20, BA3_random80)
BA3_sample <- BA3_sample[order(BA3_sample$date), ]

# BA4
set.seed(401)
BA4_first20 <- BA4_4month[1:20, ]
BA4_random80 <- BA4_4month[21:nrow(BA4_4month), ][sample(nrow(BA4_4month) - 20, 80, replace = FALSE), ]
BA4_sample <- rbind(BA4_first20, BA4_random80)
BA4_sample <- BA4_sample[order(BA4_sample$date), ]

# BA5
set.seed(502)
BA5_first20 <- BA5_4month[1:20, ]
BA5_random80 <- BA5_4month[21:nrow(BA5_4month), ][sample(nrow(BA5_4month) - 20, 80, replace = FALSE), ]
BA5_sample <- rbind(BA5_first20, BA5_random80)
BA5_sample <- BA5_sample[order(BA5_sample$date), ]

# BA6
set.seed(603)
BA6_first20 <- BA6_4month[1:20, ]
BA6_random80 <- BA6_4month[21:nrow(BA6_4month), ][sample(nrow(BA6_4month) - 20, 80, replace = FALSE), ]
BA6_sample <- rbind(BA6_first20, BA6_random80)
BA6_sample <- BA6_sample[order(BA6_sample$date), ]


## Strain name file to use for subsetting aligned fasta file
sample_list = list(BA1_sample, BA2_sample, BA3_sample, BA4_sample, BA5_sample, BA6_sample)

for (i in 1:length(sample_list)) {
  ba_sequences <- sample_list[[i]][, "seqName"] # Subset the 'seqName' column
  file_name_seq <- paste0("BA", i, "_sequences.txt") # Construct file name
  write.table(ba_sequences, file = file_name_seq, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE) # Write to file
}


## Location information
BA1_location <- BA1_sample[, c("clade_display", "pango_lineage", "seqName", "country", "province", "district")]
BA2_location <- BA2_sample[, c("clade_display", "pango_lineage", "seqName", "country", "province", "district")]
BA3_location <- BA3_sample[, c("clade_display", "pango_lineage", "seqName", "country", "province", "district")]
BA4_location <- BA4_sample[, c("clade_display", "pango_lineage", "seqName", "country", "province", "district")]
BA5_location <- BA5_sample[, c("clade_display", "pango_lineage", "seqName", "country", "province", "district")]
BA6_location <- BA6_sample[, c("clade_display", "pango_lineage", "seqName", "country", "province", "district")]

## Clean Province names and fill province_city variable
location_list = list(BA1_location, BA2_location, BA3_location, BA4_location, BA5_location, BA6_location)
#Replace correct province names (for all BA_location files)
for (i in 1:length(location_list)) {
  location_list[[i]] <- location_list[[i]] %>%
    mutate(province = case_when(
      province == "Western Cape Province" ~ "Western Cape",
      province == "North-West ZA" ~ "North West",
      province == "Kapstadt" ~ "Western Cape",
      TRUE ~ province
    ))
}

names(location_list) <- paste0("BA", 1:6, "_location") # Call the dataframe names
list2env(location_list, envir = .GlobalEnv) # Regenerate the dataframes

#Fill province_city variable for mapping (phylogeography) purposes
#Define a function to perform the transformation
mutate_location <- function(df) {
  df %>%
    mutate(district = case_when(
      is.na(district) & province == "Gauteng" ~ "Johannesburg",
      is.na(district) & province == "Western Cape" ~ "City of Cape Town",
      is.na(district) & province == "KwaZulu-Natal" ~ "eThekwini",
      is.na(district) & province == "Eastern Cape" ~ "Bhisho",
      is.na(district) & province == "Limpopo" ~ "Polokwane",
      is.na(district) & province == "Mpumalanga" ~ "Mbombela",
      is.na(district) & province == "North West" ~ "Mafikeng",
      is.na(district) & province == "Northern Cape" ~ "Kimberley",
      is.na(district) & province == "Free State" ~ "Bloemfontein",
      TRUE ~ district
    ))
}

# Apply the transformation to each dataframe in the list
for (i in 1:length(location_list)) {
  location_list[[i]] <- mutate_location(location_list[[i]])
}

# Optionally, if you want to name the list elements
names(location_list) <- paste0("BA", 1:6, "_location") # Call the dataframe names
list2env(location_list, envir = .GlobalEnv) # Regenerate the dataframes


# If geocoding fails, fix the specific location files
# BA1
BA1_location <- BA1_location %>%
  mutate(district = case_when(
    district == "",
    TRUE ~ district
  ))
# BA2
BA2_location <- BA2_location %>%
  mutate(district = case_when(
    district == "",
    TRUE ~ district
  ))
# BA3
BA3_location <- BA3_location %>%
  mutate(district = case_when(
    district == "",
    TRUE ~ district
  ))
# BA4
BA4_location <- BA4_location %>%
  mutate(district = case_when(
    district == "",
    TRUE ~ district
  ))
# BA5
BA5_location <- BA5_location %>%
  mutate(district = case_when(
    district == "",
    TRUE ~ district
  ))
# BA6
BA6_location <- BA6_location %>%
  mutate(district = case_when(
    district == "",
    TRUE ~ district
  ))


## Geocode using open street map (osmdata) library
#BA1
BA1_location <- BA1_location %>%
  mutate(
    bbox = lapply(district, getbb),
    latitude = sapply(bbox, function(b) mean(b[2,])),
    longitude = sapply(bbox, function(b) mean(b[1,]))
  )
#BA2
BA2_location <- BA2_location %>%
  mutate(
    bbox = lapply(district, getbb),
    latitude = sapply(bbox, function(b) mean(b[2,])),
    longitude = sapply(bbox, function(b) mean(b[1,]))
  )
#BA3
BA3_location <- BA3_location %>%
  mutate(
    bbox = lapply(district, getbb),
    latitude = sapply(bbox, function(b) mean(b[2,])),
    longitude = sapply(bbox, function(b) mean(b[1,]))
  )
#BA4
BA4_location <- BA4_location %>%
  mutate(
    bbox = lapply(district, getbb),
    latitude = sapply(bbox, function(b) mean(b[2,])),
    longitude = sapply(bbox, function(b) mean(b[1,]))
  )
#BA5
BA5_location <- BA5_location %>%
  mutate(
    bbox = lapply(district, getbb),
    latitude = sapply(bbox, function(b) mean(b[2,])),
    longitude = sapply(bbox, function(b) mean(b[1,]))
  )
#BA6
BA6_location <- BA6_location %>%
  mutate(
    bbox = lapply(district, getbb),
    latitude = sapply(bbox, function(b) mean(b[2,])),
    longitude = sapply(bbox, function(b) mean(b[1,]))
  )


# Plot geocoded locations on the map to check if the geocoding was performed correctly
# library(leaflet)
#BA1
leaflet(data = BA1_location) %>% addTiles() %>% addMarkers(~longitude, ~latitude)
#BA2
leaflet(data = BA2_location) %>% addTiles() %>% addMarkers(~longitude, ~latitude)
#BA3
leaflet(data = BA3_location) %>% addTiles() %>% addMarkers(~longitude, ~latitude)
#BA4
leaflet(data = BA4_location) %>% addTiles() %>% addMarkers(~longitude, ~latitude)
#BA5
leaflet(data = BA5_location) %>% addTiles() %>% addMarkers(~longitude, ~latitude)
#BA6
leaflet(data = BA6_location) %>% addTiles() %>% addMarkers(~longitude, ~latitude)


## Save location files
ba_locations = list(BA1_location, BA2_location, BA3_location, BA4_location, BA5_location, BA6_location)
for (i in 1:length(ba_locations)) {
  ba_location_final <- ba_locations[[i]][, c("seqName", "latitude", "longitude")] # Subset columns
  file_name_loc <- paste0("BA", i, "_location_file.txt") # Construct file name
  write.table(ba_location_final, file = file_name_loc, row.names = FALSE, sep = "\t", quote = FALSE) # Write to file
}


## Save full attribute information files
for (i in 1:length(ba_locations)) {
  ba_location_attribute <- ba_locations[[i]][, c("seqName", "pango_lineage", "country", "province", "district", "latitude", "longitude")]
  colnames(ba_location_attribute) <- c("strain", "pango_lineage", "country", "province", "district", "latitude", "longitude")  # Subset and rename columns
  file_name <- paste0("BA", i, "_location_attribute.txt") # Construct file name
  write.table(ba_location_attribute, file = file_name, row.names = FALSE, sep = "\t", quote = FALSE) # Write to file
}