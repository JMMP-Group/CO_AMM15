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
INLSM_DIR : str
    The path the input  LSM
INCUBE_DIR : str
    The path the input cube of bathymetry
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

import sys, platform

from datetime import datetime
import subprocess

def set_history(cube):
    cube.attributes[ 'History' ] = "Created by {} from branch {} of {} on {} ".format(script,branch[:],repos[:],current_time)
    cube.attributes[ 'Input' ]   = "GEBCO_CUBE.nc, created by  MAKE_GEBCO_CUBE.py, NWS_CUT_GEBCO_2020_TID.nc and the AMM15 unrotated coordinates file"
    cube.attributes[ 'Python version' ] = platform.python_version()
    cube.attributes[ 'System' ]  = platform.system()
    cube.attributes[ 'Release' ] = platform.release()




parser = argparse.ArgumentParser(description='Process inputs and output file paths')
parser.add_argument('-a','--AMM15_PATH', metavar='AMM15_cube_file', nargs=1,
                    help='File location of AMM15 rotated grid', required=True)
parser.add_argument('-i','--INLSM_DIR',metavar='INLSM_DIR',  nargs=1,
                    help='Path to source LSM files', required=True)
parser.add_argument('-c','--INCUBE_DIR',metavar='INCUBE_DIR',  nargs=1,
                    help='Path to source cube directory ', required=True)
parser.add_argument( '-o','--OUT_DIR',metavar='OUT_DIR', nargs=1,
                    help='Path to output to ', required=True) 

args = parser.parse_args()
if not all([args.AMM15_PATH, args.INLSM_DIR, args.INCUBE_DIR, args.OUT_DIR]):
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

print("\n     the INPUT directory for LSM as {}\n".format(  args.INLSM_DIR[0] ))
if (Path(args.INLSM_DIR[0])).is_dir():
       print(" and the {} exists.".format (args.INLSM_DIR[0]))
else:
      sys.exit("However, {} does not exist  so we exit here".format(args.INLSM_DIR[0])) 
print("\n     the INPUT directory for cube as {}\n".format(  args.INCUBE_DIR[0] ))
if (Path(args.INCUBE_DIR[0])).is_dir():
       print(" and the {} exists.".format (args.INCUBE_DIR[0]))
else:
      sys.exit("However, {} does not exist  so we exit here".format(args.INSUBE_DIR[0])) 
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


# we want to save the 2d rotated lats and lons for the extended grid in NEMO format
from get_expanded_coordinates import output_nemo_coords

output_nemo_coords(inflate_lat, inflate_lon, args.OUT_DIR[0])

# Read in the cube of gebco data this was created by MAKE_GEBCO_CUBE.py

GEBCO_RAW_cube = iris.load('{}/GEBCO_CUBE.nc'.format( args.INCUBE_DIR[0] ))[0]


# Read in the in the LSM for the GEBCO data

GEBCO_LSM_fp = Dataset('{}/NWS_CUT_GEBCO_2020_TID.nc'.format( args.INLSM_DIR[0] ),'r')
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


#=== remove land values and fill ===
GEBCO_RAW_cube.data = fill(GEBCO_RAW_cube.data)
EXTRAPOLATE_GEBCO_ON_EXPANDAMM15 = GEBCO_RAW_cube.regrid(expandcube, iris.analysis.Linear(extrapolation_mode='extrapolate'))

MASK_EXTRAPOLATE = MASK_GEBCO_ON_EXPANDAMM15


# use the masked version beyond the he AMM15 domain
# use the extrapolated version in the coee AMM15 domain

MASK_EXTRAPOLATE.data                 [extra_lat:-extra_lat,extra_lon:-extra_lon] = (
EXTRAPOLATE_GEBCO_ON_EXPANDAMM15.data [extra_lat:-extra_lat,extra_lon:-extra_lon]   )



# set global attributes to keep track of the origin of the files

#--------------------------------------------------------------------------------------
# save it to file, set up global attributes to help trace how the file was created
#--------------------------------------------------------------------------------------

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

set_history( MASK_GEBCO_ON_EXPANDAMM15 )
set_history( EXTRAPOLATE_GEBCO_ON_EXPANDAMM15 )
set_history( MASK_EXTRAPOLATE )



iris.save(MASK_GEBCO_ON_EXPANDAMM15, '{}/MASK_GEBCO_vDec2020_ON_EXPAND_AMM15.nc'.format( args.OUT_DIR[0] ))
iris.save(EXTRAPOLATE_GEBCO_ON_EXPANDAMM15, '{}/EXTRAPOLATE_GEBCO_vDec2020_ON_EXPAND_AMM15.nc'.format( args.OUT_DIR[0] ))
iris.save(MASK_EXTRAPOLATE, '{}/MASK_EXTRAPOLATE_GEBCO_vDec2020_ON_EXPAND_AMM15.nc'.format(args.OUT_DIR[0] ))





   
   
