import requests
import json

# Define the URL for the ArcGIS REST service
url = "https://services1.arcgis.com/cRvLdSPAsRupRo7I/arcgis/rest/services/US_State_Boundaries_(CONUS)/FeatureServer/0/query"

# Define the parameters for the query
params = {
    "where": "1=1",  # Query all features
    "outFields": "*",  # Get all fields
    "f": "geojson"  # Get the output in GeoJSON format
}

# Make the request to the ArcGIS REST service
response = requests.get(url, params=params)

# Check if the request was successful
if response.status_code == 200:
    # Save the response content to a GeoJSON file
    with open("US_State_Boundaries_CONUS.geojson", "w") as f:
        json.dump(response.json(), f)
    print("Data downloaded and saved as US_State_Boundaries_CONUS.geojson")
else:
    print(f"Failed to download data. HTTP status code: {response.status_code}")