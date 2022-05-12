import xarray as xr
import numpy as np
import time

import sys
import os
from dask.diagnostics import ProgressBar

import argparse
from pathlib import Path
import sys
import platform

from datetime import datetime
import subprocess

# Just parses the command line arguments in path and out path

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


infile = args.INPUT_FILE[0]

print("\n %s will be cut down to AMM15 dimensions\n"%(infile))
dsa = xr.open_mfdataset(infile,combine='by_coords').squeeze()



inflate_lat=100
inflate_lon=100


# Get core part of the domain 

print(dsa)
dscut = dsa.Bathymetry[      inflate_lat : -inflate_lat , inflate_lon : -inflate_lon ]


#dsa.keys()

LIST=list(dsa.keys())
print(LIST[0])
if(LIST[0]=='lat'):
   dlatcut = dsa.lat [      inflate_lat : -inflate_lat , inflate_lon : -inflate_lon ]
   dloncut = dsa.lon [      inflate_lat : -inflate_lat , inflate_lon : -inflate_lon ]
else:
   dlatcut = dsa.nav_lat [      inflate_lat : -inflate_lat , inflate_lon : -inflate_lon ]
   dloncut = dsa.nav_lon [      inflate_lat : -inflate_lat , inflate_lon : -inflate_lon ]




print('Storing Pure Cut Domain')

with ProgressBar():
  filename = os.path.basename(infile)
  print("filename = %s"%(filename))
  print("infile = %s"%(infile))

  dscut.to_netcdf("%s/CUTAMM15_%s"%(args.OUT_DIR[0],filename))

print(' Copying perimeter values from 4 in out to edge')

dscut.load()
for i in range(4):
    dscut[   i  ,   :  ] = dscut[   4  , :  ]
    dscut[   i  ,   :  ] = dscut[   4  , :  ]
    dscut[   :  ,   i  ] = dscut[   :  , 4  ]
    dscut[ -1-i ,   :  ] = dscut[   -5 , :  ]
    dscut[   :  , -1-i ] = dscut[   :  , -5 ]

print('Storing Cut Domain with copied outer perimieter' )

with ProgressBar():
  OUTPUT=xr.merge([dscut, dlatcut, dloncut])
  OUTPUT.attrs['description'] = 'Used Cut_and_copy_bdy_perimeter.py to cut out inner domain and apply bdy perimeter to %s'%(infile)
  OUTPUT.attrs['history'] = "Created at the Met Office, " + time.ctime(time.time())
  OUTPUT.attrs['Usage'] = "See README_AMM15_USAGE.md"
  OUTPUT.to_netcdf("%s/BDY_COPY_CUTAMM15_%s"%(args.OUT_DIR[0],filename))

print(' All Done')
