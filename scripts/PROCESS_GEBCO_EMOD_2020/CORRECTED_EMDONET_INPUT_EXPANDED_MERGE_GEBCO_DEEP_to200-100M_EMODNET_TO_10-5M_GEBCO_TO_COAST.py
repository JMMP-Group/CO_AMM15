"""
 Script to merge EMODNET and GEBCO DATA


 Uses GEBCO in the deep until 200 m
 Merge linearly between GEBCO and EMODNET from 200 m to 100m 
 Uses EDMONET from 100 m to 10 m
 Linearly interpolates EMODNET to GEBCO 10 m to 5 m
 Uses GEBCO from 5 m to the coast.
Parameters
----------
IN_DIR : str
    The path the input bathymetry
OUT_DIR : str
    The path write the merge data to

Returns
-------
Merge data
    The bathymetry merged from GEBCO and EMODNET
"""

from netCDF4 import Dataset
import numpy as np
import time

import argparse
from pathlib import Path
import sys
import platform

from datetime import datetime
import subprocess



# Just parses the command line arguments in path and out path

parser = argparse.ArgumentParser(description="Process inputs and output file paths")
parser.add_argument(
    "-e",
    "--EMOD_FILE",
    metavar="EMOD_FILE",
    nargs=1,
    help="Path to source EMODNET bathy",
    required=True,
)
parser.add_argument(
    "-g",
    "--GEB_FILE",
    metavar="GEB_FILE",
    nargs=1,
    help="Path to source GEBCO bathy",
    required=True,
)
parser.add_argument(
    "-o",
    "--OUT_DIR",
    metavar="OUT_DIR",
    nargs=1,
    help="Path to output to ",
    required=True,
)
parser.add_argument(
    "-c",
    "--COORD_FILE",
    metavar="COORD_FILE",
    nargs=1,
    help="Path to coordinates ",
    required=True,
)

args = parser.parse_args()

if not all([args.GEB_FILE, args.EMOD_FILE,args.OUT_DIR]):
    print(" Sorry All Arguments are required")
    sys.exit("Sorry all the arguments are required")


print("\n----------------------------------------------------\n")
print("Thanks, you have chosen: \n ")
print("\n     the INPUT FILE for EMOD as {}\n".format(args.EMOD_FILE[0]))
if (Path(args.EMOD_FILE[0])).is_file():
    print(" and the file {} exists.".format(args.EMOD_FILE[0]))
else:
    sys.exit("However, {} does not exist, so we exit here.".format(args.EMOD_FILE[0]))
print("\n     the INPUT FILE for GEBCO as {}\n".format(args.GEB_FILE[0]))
if (Path(args.GEB_FILE[0])).is_file():
    print(" and the file {} exists.".format(args.GEB_FILE[0]))
else:
    sys.exit("However, {} does not exist, so we exit here.".format(args.GEB_FILE[0]))
print("\n     the OUTPUT directory for processed bathy as {}\n".format(args.OUT_DIR[0]))
if (Path(args.OUT_DIR[0])).is_dir():
    print(" and the  directory {} exists.".format(args.OUT_DIR[0]))
else:
    sys.exit("However, {} does not exist, so we exit here.".format(args.OUT_DIR[0]))
print("\n   The Coordinates file selected is   {}\n".format(args.COORD_FILE[0]))
if (Path(args.COORD_FILE[0])).is_file():
    print(" and the  file {} exists.".format(args.COORD_FILE[0]))
else:
    sys.exit("However, {} does not exist, so we exit here.".format(args.COORD_FILE[0]))
print("\n----------------------------------------------------\n")

# ------------------------------------------------------------------------------
# Get input data
# ------------------------------------------------------------------------------
#%%

gebfname =('{}'.format(args.GEB_FILE[0])) 
emodfname =('{}'.format(args.EMOD_FILE[0]))
emod = Dataset(emodfname,'r')
gebc = Dataset(gebfname,'r')
emodnet_bathy = emod.variables['Bathymetry'][:,:]
gebco_bathy   = gebc.variables['Bathymetry'][:,:]
#Say 3000 to 2000

top_1 = 0200.
bot_1 = 0100.

top_2 = 0010.
bot_2 = 0005.

ramp_pass_1    = np.zeros( np.shape(emodnet_bathy) )
ramp_pass_1[:] = 1 - (emodnet_bathy[:]-bot_1 )/top_1


ramp_pass_1[ np.where( emodnet_bathy[:] >=  top_1 )] = 0.
ramp_pass_1[ np.where( emodnet_bathy[:]  <  bot_1 )] = 1.



emodnet_bath_merge_1 =  emodnet_bathy *ramp_pass_1[:] +  (1-ramp_pass_1[:])*gebco_bathy[:]


#run second pass for shallow merge
ramp_pass_2    = np.zeros(np.shape(emodnet_bath_merge_1))
ramp_pass_2[:] = 1 - (emodnet_bath_merge_1[:]-bot_2 )/top_2

ramp_pass_2[ np.where( emodnet_bath_merge_1[:] >=  top_2 )] = 0.
ramp_pass_2[ np.where( emodnet_bath_merge_1[:]  <  bot_2 )] = 1.


emodnet_bath_merge_2 =  gebco_bathy *ramp_pass_2[:] +  (1-ramp_pass_2[:])*emodnet_bath_merge_1[:]
fname2 =('{}'.format(args.COORD_FILE[0]))
coord = Dataset(fname2,'r')
lat = coord.variables['gphit'][:,:]
lon = coord.variables['glamt'][:,:]

y,x = np.shape(lat)


ncfile = Dataset('%s/CORRECTED_EXPANDED_MERGE_GEBCO_DEEP_TO_%d-%d_EMODNET_TO_%d-%d_GEBCO_TO_COAST_amm15.bathydepth.co7.cs3x.cs20.nc'%(args.OUT_DIR[0],int(top_1),int(bot_1),int(top_2),int(bot_2)),'w')
print("\n OUTPUT FILE IS:%s/CORRECTED_EXPANDED_MERGE_GEBCO_DEEP_TO_%d-%d_EMODNET_TO_%d-%d_GEBCO_TO_COAST_amm15.bathydepth.co7.cs3x.cs20.nc \n"%(args.OUT_DIR[0],int(top_1),int(bot_1),int(top_2),int(bot_2)))
y_dim = ncfile.createDimension('x', x)     # Y
x_dim = ncfile.createDimension('y', y)     # X
latout = ncfile.createVariable('lat', 'f4', ('y','x',))
lonout = ncfile.createVariable('lon', 'f4', ('y','x',))
bathy = ncfile.createVariable('Bathymetry', 'f4', ('y','x',))

latout[:]=lat
lonout[:]=lon

bathy[:]=emodnet_bath_merge_2

import os

#txt = txt.replace('\n', '')

ncfile.description = 'EXPANDED MERGE GEBCO DEEP TO %d-%d EMODNET TO %d-%d GEBCO TO COAST'%(int(top_1),int(bot_1),int(top_2),int(bot_2))

ncfile.history = "Created at the Met Office , " + time.ctime(time.time()) + " By: " +  os.path.basename(__file__) + " See README_MERGE_EXPANDED.md"

latout.units = 'degrees north'
lonout.units = 'degrees east'
bathy.units = 'meters'

ncfile.close()


