
from netCDF4 import Dataset   # Read netcdf and write netcdf
import matplotlib as mpl
mpl.use('Agg')   # If we dont want a figure to display
import numpy as np            # Get numpy
from pylab import * # gets all of matplotlib

# Need below to interpolate


import sys
#to deal with rotation load in iris
import cartopy.crs as ccrs
import iris
import iris.analysis.cartography
import numpy as np
from scipy import ndimage as nd

from mask_tools import fill

#-------------------------------------------------
# Open up the data sets
#-------------------------------------------------
base = './'

CS20_FILE = ('cs20_stat.nc')   # Read
CS20_NC = Dataset(CS20_FILE, 'r')

AMP_CS20= CS20_NC.variables['LAT'][:]
lat_CS20= CS20_NC.variables['nav_lat'][:,0]
lon_CS20= CS20_NC.variables['nav_lon'][0,:]

#For comparison read in Jenny's version of CS3X
JENNY_CS3X_FILE =('CS3X_LAT.nc')  # Read
JENNY_CS3X_NC = Dataset(JENNY_CS3X_FILE, 'r')

JENNY_AMP_CS3X= JENNY_CS3X_NC.variables['LAT'][:]
JENNY_lat_CS3X= JENNY_CS3X_NC.variables['lat'][:]
JENNY_lon_CS3X= JENNY_CS3X_NC.variables['lon'][:]



#=== remove nan values and fill ===
#Jenny version
JENNY_AMP_CS3X[JENNY_AMP_CS3X < -1.e33] = nan
JENNY_AMP_CS3X[JENNY_AMP_CS3X > 1.e33] = nan
JENNY_AMP_CS3X[JENNY_AMP_CS3X == 0] = nan
#
JENNY_AMP_CS3X = fill(JENNY_AMP_CS3X)


#CS20

AMP_CS20[AMP_CS20 == 999] = nan
AMP_CS20 = fill(AMP_CS20)

#--------------------------
# Read AMM15 unrotate
######################
# 37.5 lat , 177.5
from iris.coords import DimCoord
from iris.cube import Cube
from iris.analysis.cartography import rotate_pole

# Needs a reference AMM15 cube
AMM15_cube=iris.load('/scratch/fred/REQUIRED_INPUTS/AMM15_ROTATED_CS.nc')[0]




cs = iris.coord_systems.GeogCS(6371229)

JENNY_latitude_cs3x  = DimCoord(JENNY_lat_CS3X, standard_name='latitude', units='degrees',coord_system=cs)
JENNY_longitude_cs3x = DimCoord(JENNY_lon_CS3X, standard_name='longitude', units='degrees',coord_system=cs)

latitude_cs20  = DimCoord(lat_CS20, standard_name='latitude', units='degrees',coord_system=cs)
longitude_cs20 = DimCoord(lon_CS20, standard_name='longitude', units='degrees',coord_system=cs)


JENNY_CS3X_cube = Cube(JENNY_AMP_CS3X, standard_name='sea_floor_depth_below_geoid',units='m',
                    dim_coords_and_dims=[(JENNY_latitude_cs3x, 0), (JENNY_longitude_cs3x, 1)])

CS20_cube = Cube(AMP_CS20, standard_name='sea_floor_depth_below_geoid',units='m',
            dim_coords_and_dims=[(latitude_cs20, 0), (longitude_cs20, 1)])


rotated_CS20_amm15 = CS20_cube.regrid(AMM15_cube, iris.analysis.Linear(extrapolation_mode='mask'))
rotated_JENNY_CS3X_amm15 = JENNY_CS3X_cube.regrid(AMM15_cube, iris.analysis.Linear())
rotated_CS20_amm15 = CS20_cube.regrid(AMM15_cube, iris.analysis.Linear())




iris.save(rotated_CS20_amm15, 'COLIN_CS20_ON_AMM15.nc')
iris.save(rotated_CS20_amm15-rotated_JENNY_CS3X_amm15, 'COLIN_CS20-JENNY_CS3X_ON_AMM15.nc')


COMBINE_C3X_CS20=rotated_JENNY_CS3X_amm15.copy()
COMBINE_C3X_CS20.data[np.where(rotated_CS20_amm15.data[:]<1000.)] = rotated_CS20_amm15.data[np.where(rotated_CS20_amm15.data[:]<1000)]
COMBINE_C3X_CS20.data[np.where(np.isnan(AMM15_cube.data[:]))] = 'NaN'
iris.save(COMBINE_C3X_CS20, 'MERGE_JENNY_C3X_COLIN_CS20.nc')

