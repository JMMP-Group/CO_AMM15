# Turns out the EDMONET data as is as a repeat of lat and lon iles so we remove them here
# and output as a cube
import cartopy.crs as ccrs
import iris
import iris.analysis.cartography
import numpy as np

from scipy import ndimage as nd

#from mask_tools import fill


from iris.coords import DimCoord
from iris.cube import Cube
from iris.analysis.cartography import rotate_pole
from netCDF4 import Dataset



#EMODNET_RAW_cube=iris.load('INPUT_DATA//EMODNET_ALL_v2.nc')[0]
# Turns out lat and lon in this are not monotonic look like there is a repeat we need to clean that up

#EMODNET_RAW_fp=Dataset('INPUT_DATA//EMODNET_ALL_v2.nc','r')
EMODNET_RAW_fp=Dataset('..//MERGE/ALLmerge.nc','r')

EMODNET_LAT=EMODNET_RAW_fp.variables['lat'][:]
EMODNET_LON=EMODNET_RAW_fp.variables['lon'][:]

#EMODNET_BATHY=EMODNET_RAW_fp.variables['Bathymetry'][:]
EMODNET_BATHY=EMODNET_RAW_fp.variables['elevation'][:]
EDMONET_BATHY_B=np.copy(EMODNET_BATHY)

#Travere lats and lons to find repeats
for i in range(1,np.size(EMODNET_LON[:])):
  if(EMODNET_LON[i]<=EMODNET_LON[i-1]):
       print ('Caught a repeat at :,i,EMODNET_LON[i-1],EMODNET_LON[i],EMODNET_LON[i+1],EMODNET_LON[i+2]')
       print ('Caught a repeat at :',i,EMODNET_LON[i-1],EMODNET_LON[i],EMODNET_LON[i+1],EMODNET_LON[i+2])
EMODNET_BATHY = np.delete(EMODNET_BATHY,9482,1)
EMODNET_LON   = np.delete(EMODNET_LON,9482,0)

EMODNET_BATHY = np.delete(EMODNET_BATHY,9482,1)
EMODNET_LON   = np.delete(EMODNET_LON,9482,0)

EMODNET_BATHY = np.delete(EMODNET_BATHY,9482,1)
EMODNET_LON   = np.delete(EMODNET_LON,9482,0)


EMODNET_BATHY = np.delete(EMODNET_BATHY,18963,1)
EMODNET_LON   = np.delete(EMODNET_LON,18963,0)

EMODNET_BATHY = np.delete(EMODNET_BATHY,18963,1)
EMODNET_LON   = np.delete(EMODNET_LON,18963,0)
 
EMODNET_BATHY = np.delete(EMODNET_BATHY,18963,1)
EMODNET_LON   = np.delete(EMODNET_LON,18963,0)
 
EMODNET_BATHY = np.delete(EMODNET_BATHY,28444,1)
EMODNET_LON   = np.delete(EMODNET_LON,28444,0)
 
EMODNET_BATHY = np.delete(EMODNET_BATHY,28444,1)
EMODNET_LON   = np.delete(EMODNET_LON,28444,0)
 
EMODNET_BATHY = np.delete(EMODNET_BATHY,28444,1)
EMODNET_LON   = np.delete(EMODNET_LON,28444,0)
 
EMODNET_BATHY = np.delete(EMODNET_BATHY,28444,1)
EMODNET_LON   = np.delete(EMODNET_LON,28444,0)
 
EMODNET_BATHY = np.delete(EMODNET_BATHY,37924,1)
EMODNET_LON   = np.delete(EMODNET_LON,37924,0)
 
EMODNET_BATHY = np.delete(EMODNET_BATHY,37924,1)
EMODNET_LON   = np.delete(EMODNET_LON,37924,0)
 
EMODNET_BATHY = np.delete(EMODNET_BATHY,37924,1)
EMODNET_LON   = np.delete(EMODNET_LON,37924,0)
 
EMODNET_BATHY = np.delete(EMODNET_BATHY,47405,1)
EMODNET_LON   = np.delete(EMODNET_LON,47405,0)
 
EMODNET_BATHY = np.delete(EMODNET_BATHY,47405,1)
EMODNET_LON   = np.delete(EMODNET_LON,47405,0)
 
EMODNET_BATHY = np.delete(EMODNET_BATHY,47405,1)
EMODNET_LON   = np.delete(EMODNET_LON,47405,0)
 
EMODNET_BATHY = np.delete(EMODNET_BATHY,47405,1)
EMODNET_LON   = np.delete(EMODNET_LON,47405,0)
 
EMODNET_BATHY = np.delete(EMODNET_BATHY,56885,1)
EMODNET_LON   = np.delete(EMODNET_LON,56885,0)
 
EMODNET_BATHY = np.delete(EMODNET_BATHY,56885,1)
EMODNET_LON   = np.delete(EMODNET_LON,56885,0)
 
EMODNET_BATHY = np.delete(EMODNET_BATHY,56885,1)
EMODNET_LON   = np.delete(EMODNET_LON,56885,0)




#Repeat Check
for i in range(1,np.size(EMODNET_LON[:])):
  if(EMODNET_LON[i]<=EMODNET_LON[i-1]):
       print ('2nd check Caught a repeat at :,i,EMODNET_LON[i],EMODNET_LON[i-1],EMODNET_LON[i+1]')
       print ('2nd check Caught a repeat at :',i,EMODNET_LON[i],EMODNET_LON[i-1],EMODNET_LON[i+1])
#same for lats
for i in range(1,np.size(EMODNET_LAT[:])):
  if(EMODNET_LAT[i]<=EMODNET_LAT[i-1]):
       print ('Caught a repeat at :,i,EMODNET_LAt[i-1],EMODNET_LAt[i],EMODNET_LAT[i+1],EMODNET_LAT[i+2]')
       print ('Caught a repeat at :',i,EMODNET_LAT[i-1],EMODNET_LAT[i],EMODNET_LAT[i+1],EMODNET_LAT[i+2])
EMODNET_BATHY   = np.delete(EMODNET_BATHY,9004,0)
EMODNET_LAT   = np.delete(EMODNET_LAT,9004,0)

EMODNET_BATHY   = np.delete(EMODNET_BATHY,9004,0)
EMODNET_LAT   = np.delete(EMODNET_LAT,9004,0)

EMODNET_BATHY   = np.delete(EMODNET_BATHY,9004,0)
EMODNET_LAT   = np.delete(EMODNET_LAT,9004,0)

EMODNET_BATHY   = np.delete(EMODNET_BATHY,18005,0)
EMODNET_LAT   = np.delete(EMODNET_LAT,18005,0)

EMODNET_BATHY   = np.delete(EMODNET_BATHY,18005,0)
EMODNET_LAT   = np.delete(EMODNET_LAT,18005,0)

EMODNET_BATHY   = np.delete(EMODNET_BATHY,18005,0)
EMODNET_LAT   = np.delete(EMODNET_LAT,18005,0)

EMODNET_BATHY   = np.delete(EMODNET_BATHY,27006,0)
EMODNET_LAT   = np.delete(EMODNET_LAT,27006,0)

EMODNET_BATHY   = np.delete(EMODNET_BATHY,27006,0)
EMODNET_LAT   = np.delete(EMODNET_LAT,27006,0)

EMODNET_BATHY   = np.delete(EMODNET_BATHY,27006,0)
EMODNET_LAT   = np.delete(EMODNET_LAT,27006,0)

 
# Repeat Test
 
#same for lats
for i in range(1,np.size(EMODNET_LAT[:])):
  if(EMODNET_LAT[i]<=EMODNET_LAT[i-1]):
       print ('2nd check Caught a repeat at :,i,EMODNET_LAt[i-1],EMODNET_LAt[i],EMODNET_LAT[i+1],EMODNET_LAT[i+2]')
       print ('2nd check Caught a repeat at :',i,EMODNET_LAT[i-1],EMODNET_LAT[i],EMODNET_LAT[i+1],EMODNET_LAT[i+2])




## set up a Cube

cs = iris.coord_systems.GeogCS(6371229)
latitude_emodnet  = DimCoord(EMODNET_LAT, standard_name='latitude', units='degrees',coord_system=cs)
longitude_emodnet = DimCoord(EMODNET_LON, standard_name='longitude', units='degrees',coord_system=cs)

print (np.shape(latitude_emodnet))
print (np.shape(longitude_emodnet))
print (np.shape(EMODNET_LAT))
print (np.shape(EMODNET_LON))
print (np.shape(EMODNET_BATHY))




EMODNET_cube = Cube(EMODNET_BATHY, standard_name='sea_floor_depth_below_geoid',units='m',
            dim_coords_and_dims=[(latitude_emodnet, 0), (longitude_emodnet, 1)])


iris.save(EMODNET_cube, '../EMODNET_v2020_NO_REPEAT_LAT_LON.nc')


