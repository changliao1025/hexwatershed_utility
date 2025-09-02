import math

def to_radians(degrees):
    return degrees * math.pi / 180

def to_degrees(radians):
    return radians * 180 / math.pi

def find_great_circle_intersection(lon1, lat1, lon2, lat2, target_lon):
    """
    Find the location on the great circle that has the specified longitude.

    Args:
        lon1, lat1: Longitude and latitude of the first point (in degrees).
        lon2, lat2: Longitude and latitude of the second point (in degrees).
        target_lon: The target longitude (in degrees).

    Returns:
        (target_lon, target_lat): The longitude and latitude of the intersection point (in degrees).
    """
    # Convert coordinates to radians
    lon1, lat1, lon2, lat2, target_lon = map(to_radians, [lon1, lat1, lon2, lat2, target_lon])

    # Calculate the difference in longitudes
    d_lon = lon2 - lon1

    # Calculate the latitude of the intersection point using spherical interpolation
    Bx = math.cos(lat2) * math.cos(d_lon)
    By = math.cos(lat2) * math.sin(d_lon)
    lat_intersection = math.atan2(math.sin(lat1) + math.sin(lat2), math.sqrt((math.cos(lat1) + Bx) ** 2 + By ** 2))

    # Calculate the longitude of the intersection point
    lon_intersection = target_lon

    # Convert the intersection point back to degrees
    lon_intersection, lat_intersection = map(to_degrees, [lon_intersection, lat_intersection])

    return lon_intersection, lat_intersection

# Example usage
lon1, lat1 = -179.997342114428, 66.2769360301573
lon2, lat2 = 179.98848859037, 66.2814388694679
target_lon = 0

intersection_lon, intersection_lat = find_great_circle_intersection(lon1, lat1, lon2, lat2, target_lon)
print(f"Intersection at longitude {target_lon}: ({intersection_lon}, {intersection_lat})")