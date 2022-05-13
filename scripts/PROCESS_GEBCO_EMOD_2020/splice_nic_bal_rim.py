
import xarray as xr

import sys
import os
import subprocess
from pathlib import Path
import argparse
from datetime import datetime
import numpy as np
import platform

parser = argparse.ArgumentParser(description="Process inputs and output file paths")
parser.add_argument(
    "-i",
    "--INPUT_FILE",
    metavar="INPUT_FILE",
    nargs=1,
    help="Path to source bathy to spice on NICO bdy data",
    required=True,
)
parser.add_argument(
    "-n",
    "--NICO_FILE",
    metavar="NICO_FILE",
    nargs=1,
    help="Path to Nico Sample Bathy",
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

if (Path(args.NICO_FILE[0])).is_file():
    print(" and the file {} exists.".format(args.NICO_FILE[0]))
else:
    sys.exit("However, {} does not exist, so we exit here.".format(args.NICO_FILE[0]))

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

print("\n %s will have nico bdy baltic rim splice into it \n"%(infile))
dsa = xr.open_mfdataset(infile,combine='by_coords').squeeze()

nico = xr.open_mfdataset(args.NICO_FILE[0])

Nico  = nico.Bathymetry[:,:]
GEG  = dsa.Bathymetry[:,:]

# quick work around the differing dimension names (x,y Vs  Nlat,Nlon)
# Take a copy and then fill the data renaming  dimension names directly didn't seem to work?
Compare      = Nico.copy()
Compare.data = GEG.data

orig_minus_nico_old = Compare-Nico

orig_minus_nico_old.to_netcdf("{}/COMPARE.nc".format(args.OUT_DIR[0]))

MERGE      = GEG.copy()

MERGE.data[750:830,1447:] = Nico.data[750:830,1447:]


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


MERGE.attrs['history'] = "Created by {} from branch {} of {} on {} ".format(script,branch[:],repos[:],current_time)
MERGE.attrs['inputs']  = "{}, {}".format("TEST_SMOOTH_BDY_COPY_CUTAMM15_CORRECTED_EXPANDED_MERGE_GEBCO_DEEP_TO_200-100_EMODNET_TO_10-5_GEBCO_TO_COAST_amm15.bathydepth.co7.cs3x.cs20.nc",
                                         "NICO_MODIFIED_BATHY.nc")
MERGE.attrs['pyversion'] = platform.python_version()
MERGE.attrs['System'] = platform.system()
MERGE.attrs['Release'] = platform.release()
MERGE.attrs['Purpose'] =  "Enforce Nico Baltic Bdy Exactly"
import sys
MERGE.attrs['CommandLine'] = " ".join(sys.argv) # str(sys.argv)

filename = os.path.basename(infile)
print("Data with spliced Nico bdy baltic will be written to: {}/ADD_NICO_BALTIC_BDY_RIM_{}".format(args.OUT_DIR[0],filename))
MERGE.to_netcdf("{}/ADD_NICO_BALTIC_BDY_RIM_{}".format(args.OUT_DIR[0],filename))
