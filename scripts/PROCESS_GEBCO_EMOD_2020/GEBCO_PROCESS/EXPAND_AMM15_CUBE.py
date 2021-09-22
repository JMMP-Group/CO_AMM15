""" Interpolate GEBCO data onto Expanded AMM15 grid

We read in the rotated amm15 grid. We create a grid based on this but that is wider
so that later smoothing will not effect the bdy regions.
Here is it set arbitarliy as a 100 points extra NSWE(rotated space)

If we were just considering the AMM15 we would enforce the existing operational LSM
on this domain in the next step as well as an LAT correction derived from Surge models

However, as we go beyond the AMM15 edges we need to use a LSM as derived from GEBCO itself
We can extrapolate bathy values inside the AMM15 so that we get values in W&D areas,
when we later impose the AMM15 mask.

Thus we have a version of the data that is extrapolated into the land
and a version that has the GEBCO LSM applied directly.

Later we use the GEBCO LSM in the areas beyond the AMM15 extents and the oper LSM in the 
inner true AMM15 region


Parameters
----------
AMM15_PATH : str
    The location of the rotated AMM15 grid
IN_DIR : str
    The path the input bathymetry
OUT_DIR : str
    The path write the output cube to

Returns
-------
GEBCO interpoated on the expanded AMM15 grid

"""

import iris
from iris.cube import Cube
from iris.coords import DimCoord
from iris.analysis.cartography import rotate_pole

import numpy as np
from mask_tools import fill


from netCDF4 import Dataset


import argparse
from pathlib import Path
import sys



parser = argparse.ArgumentParser(description='Process inputs and output file paths')
parser.add_argument('-a','--AMM15_PATH', metavar='AMM15_cube_file', nargs=1,
                    help='File location of AMM15 rotated grid', required=True)
parser.add_argument('-i','--IN_DIR',metavar='IN_DIR',  nargs=1,
                    help='Path to source bathymertry files', required=True)
parser.add_argument( '-o','--OUT_DIR',metavar='OUT_DIR', nargs=1,
                    help='Path to output to ', required=True) 

args = parser.parse_args()
if not all([args.AMM15_PATH, args.IN_DIR, args.OUT_DIR]):
    print(" Sorry All Arguments are required")
    sys.exit("Sorry all the arguments are required") 


print (args.AMM15_PATH[0])
print("\n----------------------------------------------------\n")
print( "Thanks, you have chosen: \n ")
print( "      AMM15 grid file as {}\n".format (args.AMM15_PATH[0]))
if Path(args.AMM15_PATH[0]).is_file():
       print(" and the file {} exists.".format (args.AMM15_PATH[0]))
else:
      sys.exit("However, {} does not exist  so we exit here".format(args.AMM15_PATH[0])) 

print("\n     the INPUT directory for bathy as {}\n".format(  args.IN_DIR[0] ))
if (Path(args.IN_DIR[0])).is_dir():
       print(" and the {} exists.".format (args.IN_DIR[0]))
else:
      sys.exit("However, {} does not exist  so we exit here".format(args.IN_DIR[0])) 
print("\n     the OUTPUT directory for processed bathy as {}\n".format(  args.OUT_DIR[0] ))
if (Path(args.OUT_DIR[0])).is_dir():
       print(" and the {} exists.".format (args.OUT_DIR[0]))
else:
      sys.exit("However, {} does not exist  so we exit here".format(args.OUT_DIR[0])) 

#       and output for process bathy as {}.".format(args.AMM15_PATH[0], args.IN_DIR[0], args.OUT_DIR[0]) ) 
print("\n----------------------------------------------------\n")


print("\n Loading the AMM15  Cube to get base grid:\n")

AMM15_cube = iris.load(args.AMM15_PATH)[0]
X_LAT = np.array((AMM15_cube.coord('grid_latitude').points[:]))
X_LON = np.array((AMM15_cube.coord('grid_longitude').points[:]))

DEL_LAT = X_LAT[101] - X_LAT[100]
DEL_LON = X_LON[101] - X_LON[100]

extra_lat = int(100)
extra_lon = int(100)

# Start,stop, increment

start_lat = X_LAT[0]- extra_lat*DEL_LAT
start_lon = X_LON[0]- extra_lon*DEL_LON

end_lat   = start_lat + DEL_LAT * (np.size(X_LAT[:]) + 2*extra_lat)
end_lon   = start_lon + DEL_LON * (np.size(X_LON[:]) + 2*extra_lon)

inflate_lat = np.arange(start_lat,end_lat,DEL_LAT)
inflate_lon = np.arange(start_lon,end_lon,DEL_LON)

n_lat = np.size( inflate_lat )
n_lon = np.size( inflate_lon )


# Make Expansions

rotated_cs = iris.coord_systems.RotatedGeogCS(37.5, 177.5)

# set up rotated cube using AMM15 read in lats and lons

grid_latitude = DimCoord(inflate_lat,
                    standard_name='grid_latitude',
                    units='degrees',
                    coord_system=rotated_cs)

grid_longitude = DimCoord(inflate_lon,
                     standard_name='grid_longitude',
                     units='degrees',
                     coord_system=rotated_cs)

print("\n Creating the new expanded cube:\n")

expandcube = Cube(np.zeros((n_lat, n_lon), np.float32),
            dim_coords_and_dims=[(grid_latitude, 0),
                                 (grid_longitude, 1)], standard_name='sea_floor_depth_below_geoid')


print("\n Saving the expanded domain to expand_test:\n")

iris.save(expandcube, '{}/expand_test.nc'.format( args.OUT_DIR[0]) )

# Read in the cube of gebco data this was created by MAKE_GEBCO_CUBE.py

GEBCO_RAW_cube = iris.load('{}/GEBCO_CUBE.nc'.format( args.IN_DIR[0] ))[0]


# Read in the in the LSM for the GEBCO data

GEBCO_LSM_fp = Dataset('{}/NWS_CUT_GEBCO_2020_TID.nc'.format( args.IN_DIR[0] ),'r')
GEBCO_LSM = GEBCO_LSM_fp.variables['tid'][:]

#
# Fill data 
#

GEBCO_RAW_cube.data = fill(GEBCO_RAW_cube.data)

expandcube.coord_system = AMM15_cube.coord_system

GEBCO_RAW_cube.data[np.where(GEBCO_LSM.data[:] ==0) ] = np.nan

#
# masked version od GEBCO on expanded AMM15 domain
#

MASK_GEBCO_ON_EXPANDAMM15 = GEBCO_RAW_cube.regrid(expandcube, iris.analysis.Linear(extrapolation_mode='mask'))
iris.save(MASK_GEBCO_ON_EXPANDAMM15, '{}/MASK_GEBCO_vDec2020_ON_EXPAND_AMM15.nc'.format( args.OUT_DIR[0] ))


GEBCO_RAW_cube.data[np.where(GEBCO_LSM.data[:] ==0) ] = np.nan
#=== remove land values and fill ===
GEBCO_RAW_cube.data = fill(GEBCO_RAW_cube.data)
EXTRAPOLATE_GEBCO_ON_EXPANDAMM15 = GEBCO_RAW_cube.regrid(expandcube, iris.analysis.Linear(extrapolation_mode='extrapolate'))
iris.save(EXTRAPOLATE_GEBCO_ON_EXPANDAMM15, '{}/EXTRAPOLATE_GEBCO_vDec2020_ON_EXPAND_AMM15.nc'.format( args.OUT_DIR[0] ))

MASK_EXTRAPOLATE = MASK_GEBCO_ON_EXPANDAMM15

print(np.shape(MASK_EXTRAPOLATE), np.shape( EXTRAPOLATE_GEBCO_ON_EXPANDAMM15))

# use the masked version outside of the AMM15
# use the extrapolated version in the code AMM15 domain

MASK_EXTRAPOLATE.data                 [extra_lat:-extra_lat,extra_lon:-extra_lon] = (
EXTRAPOLATE_GEBCO_ON_EXPANDAMM15.data [extra_lat:-extra_lat,extra_lon:-extra_lon]   )

iris.save(MASK_EXTRAPOLATE, '{}/MASK_EXTRAPOLATE_GEBCO_vDec2020_ON_EXPAND_AMM15.nc'.format(args.OUT_DIR[0] ))

