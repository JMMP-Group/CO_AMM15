import sys, os
from   netCDF4 import Dataset, MFDataset                                # for reading/writing netcd4
#import netcdftime
import numpy as np

from scipy.interpolate import griddata

#from mpl_toolkits.basemap import Basemap, shiftgrid, interp
#from   pyroms import bathy_tools, bathy_smoothing
import bathy_tools, bathy_smoothing

import subprocess
import argparse
from pathlib import Path
from datetime import datetime
import numpy as np
import platform



#----------------------------------------------------------------------------


parser = argparse.ArgumentParser(description="Process inputs and output file paths")
parser.add_argument(
    "-i",
    "--INPUT_FILE",
    metavar="INPUT_FILE",
    nargs=1,
    help="Path to source bathy to be cut",
    required=True,
)
parser.add_argument(
    "-o",
    "--OUT_DIR",
    metavar="OUT_DIR",
    nargs=1,
    help="Path to output",
    required=True,
)

args = parser.parse_args()


print("\n----------------------------------------------------\n")
print("Thanks, you have chosen: \n ")
print("\n     the INPUT FILE as {}\n".format(args.INPUT_FILE[0]))
if (Path(args.INPUT_FILE[0])).is_file():
    print(" and the file {} exists.".format(args.INPUT_FILE[0]))
else:
    sys.exit("However, {} does not exist, so we exit here.".format(args.INPUT_FILE[0]))

print("\n     the OUTPUT directory for processed bathy as {}\n".format(args.OUT_DIR[0]))
if (Path(args.OUT_DIR[0])).is_dir():
    print(" and the  directory {} exists.".format(args.OUT_DIR[0]))
else:
    sys.exit("However, {} does not exist, so we exit here.".format(args.OUT_DIR[0]))

# ------------------------------------------------------------------------------
# Get input data
# ------------------------------------------------------------------------------
#%%


srcFile = args.INPUT_FILE[0]
filename = os.path.basename(srcFile)

outFile= "%s/SMOOTH_%s"%(args.OUT_DIR[0],filename)


## READ SRC BATHYMETRY
nc_src  = Dataset( srcFile, 'r' )
lon_src = nc_src.variables[ 'lon' ][:]
lat_src = nc_src.variables[ 'lat' ][:]
bat_src = nc_src.variables[ 'Bathymetry' ][:]

#EOD deal with nans
bat_src[np.where(~np.isfinite(bat_src))] = -11
print (srcFile, "loaded", lon_src.shape, bat_src.shape)
bat_src[:]=bat_src+10 # Makes ponts we consider ok
indbat = ( bat_src <= 0.)
## EOD this is not true for WAD say set some limit like 10

## NETCDF OUTPUT
ncout = Dataset( outFile, 'w', format='NETCDF3_CLASSIC' )
ncout.createDimension( 'lat', bat_src.shape[0] )
ncout.createDimension( 'lon', bat_src.shape[1] )
nlon = ncout.createVariable( 'nlon', 'f4', ('lat', 'lon',), zlib='True' )
nlat = ncout.createVariable( 'nlat', 'f4', ('lat', 'lon',), zlib='True' )
nlon[:]  = lon_src; nlat[:] = lat_src
bout    = ncout.createVariable( "Maked_Bathymetry_old", 'f4', ('lat','lon',), zlib='True', fill_value=-999. )
cout    = ncout.createVariable( "Bathymetry_old", 'f4', ('lat','lon',), zlib='True', fill_value=-999. )

print ("SRC", np.min(bat_src), np.max(bat_src))

## READ COORDINATES AND EXTRACT FIRST EDGE
mask = np.ones_like( bat_src  )
maskb = np.zeros_like( bat_src  )
#ind = ( bat_src <= 0. ); mask[ind] = 0  ## place land
ind = ( bat_src == -1 ); mask[ind] = 0  ## place land
# Set a Lon limit and Lat limit
#maskb[800:1050,1400:] = 1
maskb[700:950,1300:] = 1

mask=mask*maskb


#masked_bat_src = np.ma.array ( bat_src, mask=(mask==-11.))
masked_bat_src = np.ma.array ( bat_src, mask=(bat_src==-1))
masked_bat_src = np.ma.array ( masked_bat_src, mask=(maskb==0))
bout[:] = masked_bat_src-10
cout[:] = bat_src-10
#EODind = ( bat_src <= -10. ); mask[ind] = 0  ## place land
#ind = ( bat_src >= 75.); mask[ind] = 0  ## place land over deep water to speed up
#Do the whole thing now we dont have issues with precision
#ind = ( bat_src >= 10+75.); mask[ind] = 0  ## place land over deep water to speed up


print ("mask created", mask.shape, np.min(mask), np.max(mask))
mout    = ncout.createVariable( "LSMMask", 'f4', ('lat','lon',), zlib='True', fill_value=-999. )
mout[:] = mask
moutb    = ncout.createVariable( "BALTICMask", 'f4', ('lat','lon',), zlib='True', fill_value=-999. )
#
moutb[:] = maskb*(bat_src)

# Check bathymetry roughness
RoughMat = bathy_tools.RoughnessMatrix( bat_src, mask)
print ('Ini Max Roughness Value: {0:.3f}'.format( RoughMat.max() ))

# Smooth the raw bathy using the direct iterative method from Martinho and Batteen (2006)
rx0_max = 0.1
if rx0_max != None : 
   bat_new = bathy_smoothing.smoothing_Positive_rx0( mask, bat_src, rx0_max )
   # check bathymetry roughness again
   RoughMat = bathy_tools.RoughnessMatrix( bat_new, mask )
   print ('New Max Roughness Value: {0:.3f}'.format( RoughMat.max() ))
bat_new=bat_new-10.0
bat_src=bat_src-10.0

bat_new[indbat] = -11.
bat_src[indbat] = -11.


bout2    = ncout.createVariable( "Bathymetry", 'f4', ('lat','lon',), zlib='True', fill_value=-999. )
masked_bat_new = np.ma.array ( bat_new, mask=(bat_new==-11))
bout2[:] = masked_bat_new

anom    =  ncout.createVariable( "Anomaly", 'f4', ('lat','lon',), zlib='True', fill_value=-999. )
anom[:] = bat_new - bat_src


now = datetime.now()
current_time = now.strftime("%Y/%m/%d %H:%M:%S")
repos = subprocess.run(
    ["git", "config", "--get", "remote.origin.url"],
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    check=True
)
repos = repos.stdout.decode("utf-8").strip("\n")
branch = subprocess.run(
    ["git", "rev-parse", "--abbrev-ref", "HEAD"],
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    check=True
)
branch = branch.stdout.decode("utf-8").strip("\n")
parser = argparse.ArgumentParser(description="Process inputs and output file paths")
script = parser.prog



ncout.history = "Created by {} from branch {} of {} on {} ".format(script,branch[:],repos[:],current_time)
ncout.inputs  = "{}, {}".format(srcFile,outFile)
ncout.pyversion = platform.python_version()
ncout.System = platform.system()
ncout.Release = platform.release()
ncout.Purpose =  "Smooth to rmax = {} in Baltic area 7[00:950,1300:]".format(rx0_max)
import sys
ncout.CommandLine = " ".join(sys.argv) # str(sys.argv)

ncout.close()

bat_src = np.ma.array ( bat_src, mask=(mask==-11.))
bat_new = np.ma.array ( bat_new, mask=(mask==-11.))














