""" Convert Raw bathy into Cube Format

We use the iris interp to later map the input bathy to AMM15.
To do that we convert the bathy into cube format here

Parameters
----------
IN_DIR : str
    The path the input bathymetry
OUT_DIR : str
    The path write the output cube to

Returns
-------
iris cube
    The bathymetry in cube format
"""

import iris
import numpy as np

from iris.coords import DimCoord
from iris.cube import Cube
from netCDF4 import Dataset

import argparse
from pathlib import Path
import sys, platform

from datetime import datetime
import subprocess




# Just parses the command line arguments in path and out path

parser = argparse.ArgumentParser(description='Process inputs and output file paths')
parser.add_argument('-i','--IN_DIR',metavar='IN_DIR',  nargs=1,
                    help='Path to source bathymertry files', required=True)
parser.add_argument( '-o','--OUT_DIR',metavar='OUT_DIR', nargs=1,
                    help='Path to output to ', required=True) 

args = parser.parse_args()

if not all([args.IN_DIR, args.OUT_DIR]):
    print(" Sorry All Arguments are required")
    sys.exit("Sorry all the arguments are required") 


print("\n----------------------------------------------------\n")
print( "Thanks, you have chosen: \n ")
print("\n     the INPUT directory for bathy as {}\n".format(  args.IN_DIR[0] ))
if (Path(args.IN_DIR[0])).is_dir():
       print(" and the directory{} exists.".format (args.IN_DIR[0]))
else:
      sys.exit("However, {} does not exist, so we exit here.".format(args.IN_DIR[0])) 
print("\n     the OUTPUT directory for processed bathy as {}\n".format(  args.OUT_DIR[0] ))
if (Path(args.OUT_DIR[0])).is_dir():
       print(" and the  directory {} exists.".format (args.OUT_DIR[0]))
else:
      sys.exit("However, {} does not exist, so we exit here.".format(args.OUT_DIR[0])) 
print("\n----------------------------------------------------\n")

#------------------------------------------------------------------------------
# Get input data
#------------------------------------------------------------------------------

GEBCO_RAW_fp = Dataset('{}/NWS_CUT_GEBCO.nc'.format(args.IN_DIR[0]),'r')

GEBCO_LAT = GEBCO_RAW_fp.variables['lat'][:]
GEBCO_LON = GEBCO_RAW_fp.variables['lon'][:]

GEBCO_BATHY = GEBCO_RAW_fp.variables['elevation'][:]

GEBCO_BATHY = GEBCO_BATHY.astype(float)

EDMONET_BATHY_B = np.copy(GEBCO_BATHY)

#--------------------------------------------------------------------------------------
# set up the cube
#--------------------------------------------------------------------------------------

cs = iris.coord_systems.GeogCS(6371229)

latitude_emodnet  = DimCoord(GEBCO_LAT, standard_name='latitude', units='degrees',coord_system=cs)
longitude_emodnet = DimCoord(GEBCO_LON, standard_name='longitude', units='degrees',coord_system=cs)


GEBCO_cube = Cube(GEBCO_BATHY, standard_name='sea_floor_depth_below_geoid',units='m',
            dim_coords_and_dims=[(latitude_emodnet, 0), (longitude_emodnet, 1)])

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

GEBCO_cube.attributes[ 'History' ] = "Created by {} from branch {} of {} on {} ".format(script,branch[:],repos[:],current_time)
GEBCO_cube.attributes[ 'Input' ] = "NWS_CUT_GEBCO.nc"
GEBCO_cube.attributes[ 'Python version' ] = platform.python_version()
GEBCO_cube.attributes[ 'System' ] = platform.system()
GEBCO_cube.attributes[ 'Release' ] = platform.release()


iris.save(GEBCO_cube, '{}/GEBCO_CUBE.nc'.format(args.OUT_DIR[0]))


