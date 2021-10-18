""" Convert Raw bathy into Cube Format

We use the iris interp to later map the input bathy to AMM15.
To do that we convert the bathy into cube format here

Note mergin the tiled  EDMONET data reults in a  repeat of lat and lon lines 
so we remove them here and output as a cube

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
import xarray as xr



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

#EMODNET_RAW_fp = Dataset('{}/ALLmerge.nc'.format(args.IN_DIR[0]),'r')

#EMODNET_LAT = EMODNET_RAW_fp.variables['lat'][:]
#EMODNET_LON = EMODNET_RAW_fp.variables['lon'][:]

#EMODNET_BATHY = EMODNET_RAW_fp.variables['elevation'][:]

EMODNET_RAW = xr.open_mfdataset( '{}/ALLmerge.nc'.format(args.IN_DIR[0]), parallel=True)
EMODNET_LAT = EMODNET_RAW.lat[:]
EMODNET_LON = EMODNET_RAW.lon[:]
EMODNET_BATHY = EMODNET_RAW.elevation[:]

print(EMODNET_LAT)
print(EMODNET_BATHY)
#EDMONET_BATHY_B = np.copy(EMODNET_BATHY)

print ("Data Read In")

#Travere lats and lons to find repeats
for i in range(1,np.size(EMODNET_LON[:])):
  if(EMODNET_LON[i] <= EMODNET_LON[i-1]):
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


print ("Before Repeat Check")


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
print ("set up a Cube")

cs = iris.coord_systems.GeogCS(6371229)
latitude_emodnet  = DimCoord(EMODNET_LAT, standard_name='latitude', units='degrees',coord_system=cs)
longitude_emodnet = DimCoord(EMODNET_LON, standard_name='longitude', units='degrees',coord_system=cs)



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



EMODNET_cube = Cube(EMODNET_BATHY, standard_name='sea_floor_depth_below_geoid',units='m',
            dim_coords_and_dims=[(latitude_emodnet, 0), (longitude_emodnet, 1)])
EMODNET_cube.attributes[ 'History' ] = "Created by {} from branch {} of {} on {} ".format(script,branch[:],repos[:],current_time)
EMODNET_cube.attributes[ 'Input' ] = "Allmerge.nc"
EMODNET_cube.attributes[ 'Python version' ] = platform.python_version()
EMODNET_cube.attributes[ 'System' ] = platform.system()
EMODNET_cube.attributes[ 'Release' ] = platform.release()


iris.save(EMODNET_cube, '{}/EMODNET_v2020_NO_REPEAT_LAT_LON.nc'.format(args.OUT_DIR[0]))


