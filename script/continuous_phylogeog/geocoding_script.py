import pandas as pd
import googlemaps
import time

# Load the data
location_metadata_path = 'location_metadata.tsv'
most_populated_cities_path = 'South_Africa_Provinces_Most_Populated_Cities.csv'

location_metadata = pd.read_csv(location_metadata_path, sep='\t')
most_populated_cities = pd.read_csv(most_populated_cities_path)

# Split into rows with and without a location specified
with_location = location_metadata[location_metadata['location'].notna()].copy()
without_location = location_metadata[location_metadata['location'].isna()].copy()

# Merge without_location with most populated cities data on division/province
without_location = without_location.merge(
    most_populated_cities,
    left_on='division',
    right_on='Province',
    how='left'
)
without_location = without_location[['strain', 'country', 'division', 'City', 'Latitude', 'Longitude']]
without_location.rename(columns={'City': 'location'}, inplace=True)

# Initialize Google Maps client
gmaps = googlemaps.Client(key='')

# Function to geocode a location
def geocode_location(location):
    try:
        geocode_result = gmaps.geocode(location)
        if geocode_result:
            location_data = geocode_result[0]['geometry']['location']
            return location_data['lat'], location_data['lng']
        return None, None
    except Exception as e:
        print(f"Error geocoding {location}: {e}")
        return None, None

# Apply geocoding to rows with location
with_location['Latitude'], with_location['Longitude'] = zip(
    *with_location['location'].apply(geocode_location)
)

# Combine geocoded data with fallback
final_data = pd.concat([with_location, without_location], ignore_index=True)

# Save the final data to a CSV file
final_data.to_csv('locations_metadata_geocoded.tsv', sep='\t', index=False)

# Extract only the strain, latitude, and longitude columns
location_file = final_data[['strain', 'Latitude', 'Longitude']].copy()

# Rename columns for lowercase latitude and longitude
location_file.rename(columns={'Latitude': 'latitude', 'Longitude': 'longitude'}, inplace=True)

# Save the strain, latitude, and longitude to a separate TSV file
location_file.to_csv('location_file.tsv', sep='\t', index=False)

print("Geocoding complete. Results saved to 'locations_metadata_geocoded.tsv'.")
print("Strain, latitude, and longitude saved to 'location_file.tsv'.")