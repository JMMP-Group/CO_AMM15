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






import argparse
from pathlib import Path
import sys
import platform

from datetime import datetime
import subprocess

import dask
from dask.diagnostics import ProgressBar

import iris
from iris.coords import DimCoord
from iris.cube import Cube
import numpy as np
import xarray as xr
#%%

# Just parses the command line arguments in path and out path

parser = argparse.ArgumentParser(description="Process inputs and output file paths")
parser.add_argument(
    "-i",
    "--IN_DIR",
    metavar="IN_DIR",
    nargs=1,
    help="Path to source bathymertry files",
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

args = parser.parse_args()

if not all([args.IN_DIR, args.OUT_DIR]):
    print(" Sorry All Arguments are required")
    sys.exit("Sorry all the arguments are required")


print("\n----------------------------------------------------\n")
print("Thanks, you have chosen: \n ")
print("\n     the INPUT directory for bathy as {}\n".format(args.IN_DIR[0]))
if (Path(args.IN_DIR[0])).is_dir():
    print(" and the directory{} exists.".format(args.IN_DIR[0]))
else:
    sys.exit("However, {} does not exist, so we exit here.".format(args.IN_DIR[0]))
print("\n     the OUTPUT directory for processed bathy as {}\n".format(args.OUT_DIR[0]))
if (Path(args.OUT_DIR[0])).is_dir():
    print(" and the  directory {} exists.".format(args.OUT_DIR[0]))
else:
    sys.exit("However, {} does not exist, so we exit here.".format(args.OUT_DIR[0]))
print("\n----------------------------------------------------\n")

# ------------------------------------------------------------------------------
# Get input data
# ------------------------------------------------------------------------------
#%%
EMODNET_RAWA = xr.open_mfdataset("{}/ALLmerge.nc".format(args.IN_DIR[0]), parallel=True)
EMODNET_LAT = EMODNET_RAWA.lat[:]
EMODNET_LON = EMODNET_RAWA.lon[:]

EMODNET_RAWB = xr.open_mfdataset(
    "{}/ALLmerge.nc".format(args.IN_DIR[0]),
    parallel=True,
    chunks=({"lat": 100, "lon": -1}),
)
EMODNET_BATHY = EMODNET_RAWB.elevation[:]

print(EMODNET_LAT)
print(EMODNET_BATHY)
# EDMONET_BATHY_B = np.copy(EMODNET_BATHY)

print("Data Read In")


#%%
ind_x = np.zeros(1, dtype=int)
i = 1
while i < np.size(EMODNET_LON.data[:]):
    if (EMODNET_LON.data[i] - EMODNET_LON.data[i - 1]) < 1e-10:
        k = i - 1
        print("caught a repeat", EMODNET_LON.data[i], EMODNET_LON.data[i - 1])
        while (EMODNET_LON.data[i] - EMODNET_LON.data[k]) < 1e-10:
            print(
                "Caught a further repeat",
                EMODNET_LON[i].data,
                EMODNET_LON.data[k],
                -EMODNET_LON.data[i],
                EMODNET_LON.data[k],
            )
            i = i + 1
    ind_x = np.append(ind_x[:], i)
    i = i + 1
print(ind_x)

ind_x = xr.DataArray(ind_x, dims=["lon"])
#%%
ind_y = np.zeros(1, dtype=int)
i = 1
while i < np.size(EMODNET_LAT.data[:]):
    if (EMODNET_LAT.data[i] - EMODNET_LAT.data[i - 1]) < 1e-10:
        k = i - 1
        print("caught a repeat", EMODNET_LAT.data[i], EMODNET_LAT.data[i - 1])
        while (EMODNET_LAT.data[i] - EMODNET_LAT.data[k]) < 1e-10:
            print(
                "Caught a further repeat",
                EMODNET_LAT[i].data,
                EMODNET_LAT.data[k],
                -EMODNET_LAT.data[i],
                EMODNET_LAT.data[k],
            )
            i = i + 1
    ind_y = np.append(ind_y[:], i)
    i = i + 1
print(ind_y)
ind_y = xr.DataArray(ind_y, dims=["lat"])

#%%

dask.config.set(**{"array.slicing.split_large_chunks": False})
EMODNET_BATHY = EMODNET_BATHY[ind_y, ind_x]


## set up a Cube
print("set up a Cube")

# Rename to eb compatible with expanded cube later
EMODNET_BATHY = EMODNET_BATHY.rename("sea_floor_depth_below_geoid")

#%%
#
#%%
now = datetime.now()
current_time = now.strftime("%Y/%M/%d %H:%M:%S")
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
script = parser.prog
#%%
cs = iris.coord_systems.GeogCS(6371229)
latitude_emodnet = DimCoord(
    EMODNET_LAT[ind_y], standard_name="latitude", units="degrees", coord_system=cs
)
longitude_emodnet = DimCoord(
    EMODNET_LON[ind_x], standard_name="longitude", units="degrees", coord_system=cs
)
#%%

#%%
EMODNET_cube = EMODNET_BATHY.to_iris()
EMODNET_cube.standard_name = "sea_floor_depth_below_geoid"
#%%
EMODNETB_cube = Cube(
    EMODNET_BATHY,
    standard_name="sea_floor_depth_below_geoid",
    units="m",
    dim_coords_and_dims=[(latitude_emodnet, 0), (longitude_emodnet, 1)],
)

EMODNET_cube.coord_system = EMODNETB_cube.coord_system
#%%

EMODNET_cube.attributes["History"] = "Created by {} from branch {} of {} on {} ".format(
    script, branch[:], repos[:], current_time
)
EMODNET_cube.attributes["Input"] = "Allmerge.nc"
EMODNET_cube.attributes["Python version"] = platform.python_version()  #
EMODNET_cube.attributes["System"] = platform.system()
EMODNET_cube.attributes["Release"] = platform.release()
#%%
with ProgressBar():
    iris.save(
        EMODNET_cube, "{}/EMODNET_v2020_NO_REPEAT_LAT_LON.nc".format(args.OUT_DIR[0])
    )
