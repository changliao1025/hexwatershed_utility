import math
import numpy as np
def xyz_to_lonlat(x, y, z):
    # Normalize the coordinates
    norm = np.sqrt(x**2 + y**2 + z**2)
    x /= norm
    y /= norm
    z /= norm

    # Convert to spherical coordinates
    lon = np.arctan2(y, x)
    lat = np.arcsin(z)

    # Convert radians to degrees
    lon = np.degrees(lon)
    lat = np.degrees(lat)

    return float(lon), float(lat)
class Pnt:
    def __init__(self, x=0.0, y=0.0, z=0.0, idx=0, positiveLonRange=True):
        self.x = x
        self.y = y
        self.z = z
        self.idx = idx
        self.positiveLonRange = positiveLonRange
        self.buildLat()
        self.buildLon()

    def normalize(self):
        norm = self.x * self.x + self.y * self.y + self.z * self.z
        if norm == 0:
            print("Pnt: normalize")
            print(self.x, self.y, self.z, self.idx)
            assert norm != 0, "Normalization error: norm is zero"
        norm = math.sqrt(norm)

        self.x = self.x / norm
        self.y = self.y / norm
        self.z = self.z / norm

    def buildLat(self):
        dl = math.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)
        self.lat = math.asin(self.z / dl)

    def buildLon(self):
        lon = math.atan2(self.y, self.x)

        # If the prime meridian is the minimum, this means the degree
        # range is 0-360, so we need to translate.
        if self.positiveLonRange:
            if lon < 0.0:
                lon = 2.0 * math.pi + lon

        self.lon = lon

#do some test,
#create a point at the north pole
p = Pnt(0, 0, 1, 0)

#normalize the point
p.normalize()
#get lon an  lat
dLongitude_r =  float(p.lon)
dLatitude_r =  float(p.lat)
dLongitude360 =  dLongitude_r / math.pi * 180
dLongitude = dLongitude360
dLatitude =  dLatitude_r / math.pi * 180
print(dLongitude, dLatitude)

#convert the point to lon lat using the xyz_to_lonlat function
dLongitude, dLatitude = xyz_to_lonlat(p.x, p.y, p.z)
print(dLongitude, dLatitude)

#use a different cite location such as Chicago
#set x y z to the location of Chicago on the earth
#set chicago lon and lat
#Chicago, IL, USA
#41.8781° N, 87.6298° W
lon = -87.6298
lat = 41.8781
#convert lon lat to xyz
x = math.cos(math.radians(lat)) * math.cos(math.radians(lon))
y = math.cos(math.radians(lat)) * math.sin(math.radians(lon))
z = math.sin(math.radians(lat))
#use earth radius to scale the x y z
earth_radius = 6371.0
x *= earth_radius
y *= earth_radius
z *= earth_radius


dLongitude, dLatitude = xyz_to_lonlat(x, y, z)
print(dLongitude, dLatitude)

p = Pnt(x, y, z, 0)

#normalize the point
p.normalize()
#get lon an  lat
dLongitude_r =  float(p.lon)
dLatitude_r =  float(p.lat)
dLongitude360 =  dLongitude_r / math.pi * 180
dLongitude = dLongitude360
dLatitude =  dLatitude_r / math.pi * 180
print(dLongitude, dLatitude)
