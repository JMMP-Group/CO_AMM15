""" Combine CS20, CS3X LAT data into one data set for AMM15 grid.



We read in the two surgae dat set, CS20 , CS3X interpolate each to the AMM15 grid 
and combineinto one data set. CS3X goes out to the old AMM7 bdy limits CS20 only as far the shelf break
Parameters
----------
IN_DIR : str
    The path to the surge datasets and AMM15 grid
OUT_DIR : str
    The path write the output cube to

Returns
-------
LAT correction as a combination of CS3X and CS20 dat aon the AMM15 grid

"""



from netCDF4 import Dataset   # Read netcdf and write netcdf
import numpy as np            # Get numpy

# Need below to interpolate


#to deal with rotation load in iris
import cartopy.crs as ccrs
import iris
import iris.analysis.cartography


import argparse
from pathlib import Path

import sys, platform
sys.path.append(r'../')
from mask_tools import fill


from datetime import datetime
import subprocess

def set_history(cube):
    cube.attributes[ 'History' ] = "Created by {} from branch {} of {} on {} ".format(script,branch[:],repos[:],current_time)
    cube.attributes[ 'Input' ]   = "{}".format(args.IN_DIR[0])
    cube.attributes[ 'Python version' ] = platform.python_version()
    cube.attributes[ 'System' ]  = platform.system()
    cube.attributes[ 'Release' ] = platform.release()



# Parse arguments

parser = argparse.ArgumentParser(description='Process inputs and output file paths')
parser.add_argument('-i','--IN_DIR',metavar='IN_DIR',  nargs=1,
                    help='Path to input files', required=True)
parser.add_argument( '-o','--OUT_DIR',metavar='OUT_DIR', nargs=1,
                    help='Path to output to ', required=True)

args = parser.parse_args()
if not all([ args.IN_DIR, args.OUT_DIR]):
    print(" Sorry All Arguments are required")
    sys.exit("Sorry all the arguments are required")


print("\n----------------------------------------------------\n")
print( "Thanks, you have chosen: \n ")

print("\n     the INPUT directory as {}\n".format(  args.IN_DIR[0] ))
if (Path(args.IN_DIR[0])).is_dir():
       print(" and the dir{} exists.".format (args.IN_DIR[0]))
else:
      sys.exit("However, {} does not exist, so we exit here".format(args.IN_DIR[0]))
print("\n     the OUTPUT directory for processed bathy as {}\n".format(  args.OUT_DIR[0] ))
if (Path(args.OUT_DIR[0])).is_dir():
       print(" and the {} exists.".format (args.OUT_DIR[0]))
else:
      sys.exit("However, {} does not exist  so we exit here".format(args.OUT_DIR[0]))

print("\n----------------------------------------------------\n")




#-------------------------------------------------
# Open up the data sets
#-------------------------------------------------


CS20_FILE = ('{}/cs20_stat.nc'.format(args.IN_DIR[0]))   # Read
CS20_NC = Dataset(CS20_FILE, 'r')

AMP_CS20= CS20_NC.variables['LAT'][:]
lat_CS20= CS20_NC.variables['nav_lat'][:,0]
lon_CS20= CS20_NC.variables['nav_lon'][0,:]

#For comparison read in Jenny's version of CS3X
CS3X_FILE =('{}/CS3X_LAT.nc'.format(args.IN_DIR[0]))  # Read
CS3X_NC = Dataset(CS3X_FILE, 'r')

AMP_CS3X= CS3X_NC.variables['LAT'][:]
lat_CS3X= CS3X_NC.variables['lat'][:]
lon_CS3X= CS3X_NC.variables['lon'][:]



#=== remove nan values and fill ===
AMP_CS3X[AMP_CS3X < -1.e33] = np.nan
AMP_CS3X[AMP_CS3X > 1.e33] = np.nan
AMP_CS3X[AMP_CS3X == 0] = np.nan
#
AMP_CS3X = fill(AMP_CS3X)


#CS20

AMP_CS20[AMP_CS20 == 999] = np.nan
AMP_CS20[AMP_CS20 < -999] = np.nan
AMP_CS20 = fill(AMP_CS20)

#--------------------------
# Read AMM15 unrotate
######################
# 37.5 lat , 177.5
from iris.coords import DimCoord
from iris.cube import Cube
from iris.analysis.cartography import rotate_pole

# Needs a reference AMM15 cube
AMM15_cube=iris.load('{}/AMM15_ROTATED_CS.nc'.format(args.IN_DIR[0]))[0]




cs = iris.coord_systems.GeogCS(6371229)

latitude_cs3x  = DimCoord(lat_CS3X, standard_name='latitude', units='degrees',coord_system=cs)
longitude_cs3x = DimCoord(lon_CS3X, standard_name='longitude', units='degrees',coord_system=cs)

latitude_cs20  = DimCoord(lat_CS20, standard_name='latitude', units='degrees',coord_system=cs)
longitude_cs20 = DimCoord(lon_CS20, standard_name='longitude', units='degrees',coord_system=cs)


CS3X_cube = Cube(AMP_CS3X, standard_name='sea_floor_depth_below_geoid',units='m',
                    dim_coords_and_dims=[(latitude_cs3x, 0), (longitude_cs3x, 1)])

CS20_cube = Cube(AMP_CS20, standard_name='sea_floor_depth_below_geoid',units='m',
            dim_coords_and_dims=[(latitude_cs20, 0), (longitude_cs20, 1)])


rotated_CS20_amm15 = CS20_cube.regrid(AMM15_cube, iris.analysis.Linear(extrapolation_mode='mask'))
rotated_CS3X_amm15 = CS3X_cube.regrid(AMM15_cube, iris.analysis.Linear())
rotated_CS20_amm15 = CS20_cube.regrid(AMM15_cube, iris.analysis.Linear())



now = datetime.now()
current_time = now.strftime("%Y/%M/%d %H:%M:%S")

repos = subprocess.run(['git', 'config', '--get', 'remote.origin.url'],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
repos = repos.stdout.decode('utf-8').strip('\n')

branch = subprocess.run(['git', 'rev-parse', '--abbrev-ref',  'HEAD'],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
branch = branch.stdout.decode('utf-8').strip('\n')

script = parser.prog






set_history( rotated_CS20_amm15 )
iris.save(rotated_CS20_amm15, '{}/COLIN_CS20_ON_AMM15.nc'.format(args.OUT_DIR[0]))
X=rotated_CS20_amm15-rotated_CS3X_amm15
set_history( X )
iris.save(X, '{}/COLIN_CS20-CS3X_ON_AMM15.nc'.format(args.OUT_DIR[0]))



COMBINE_C3X_CS20 = rotated_CS3X_amm15.copy()
COMBINE_C3X_CS20.data[np.where(rotated_CS20_amm15.data[:]<1000.)] = rotated_CS20_amm15.data[np.where(rotated_CS20_amm15.data[:]<1000)]
COMBINE_C3X_CS20.data[np.where(np.isnan(AMM15_cube.data[:]))] = 'NaN'

set_history( COMBINE_C3X_CS20 )
iris.save(COMBINE_C3X_CS20, '{}/MERGE_CS3X_COLIN_CS20.nc'.format(args.OUT_DIR[0]))

