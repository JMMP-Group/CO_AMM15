""" Interpolate GEBCO data onto Expanded AMM15 grid
We read in the rotated amm15 grid. We create a grid based on this but that is wider
so that later smoothing will not effect the bdy regions.
Here is it set arbitarliy as a 100 points extra NSWE (rotated space)

If we were just considering the AMM15 we would enforce the existing operational LSM
on this domain in the next step as well as an LAT correction derived from Surge models

However, as we go beyond the AMM15 edges we need to use a LSM as derived from GEBCO itself
We can extrapolate bathy values inside the AMM15 so that we get values in W&D areas,
when we later impose the AMM15 mask.

Thus we have a version of the data that is extrapolated into the land
and a version that has the EMODET LSM applied directly.

Later we use the EMODET LSM in the areas beyond the AMM15 extents and the oper LSM in the
inner true AMM15 region

Note on Spice we can use enough RAM to process the iris regrid all in one but on
more reasonable machines is is not possible and iris regrid cannot chunk in the horizontal

So we have a work aroudn where we make a small section of destination grid
and a small section of src data that covers the destination gris plus extra (1000)
so the infilling is done consistently

then at the end we can rejoin all the interpolated subsections into a single file


Parameters
----------
AMM15_PATH : str
    The location of the rotated AMM15 grid
INCUBE_DIR : str
    The path the input cube of bathymetry
INLSM_DIR : str
    The path the input  LSM
OUT_DIR : str
    The path write the output cube to

Returns
-------
GEBCO interpolated on the expanded AMM15 grid

"""
#%%


import os
import glob

from netCDF4 import Dataset


import argparse
from pathlib import Path

import sys
import platform
from datetime import datetime
import subprocess


import iris
import numpy as np
from iris.cube import Cube
from iris.coords import DimCoord
from dask.diagnostics import ProgressBar
import xarray as xr


sys.path.append(r"../")
from mask_tools import fill
from get_expanded_coordinates import output_nemo_coords


def find_nearest(array, value):
    """
    Finds the nearest value and index of that value in an ordered array to a specifed value

             Parameters:
                     array[] (float): ordered array of values (lats/lons)
                     value (float): Desired point

             Returns:
                     array[idx] (float): value at array nearest to value
                     idx(int) : index of array closes to value

    """

    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx


def set_history(cube):
    """
    Adds a record to the cubes history of when and how it was made

            Parameters:
                    cube(iris cube) : cube to be stored in netcdf
    """
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

    script = parser.prog
    cube.attributes["History"] = "Created by {} from branch {} of {} on {} ".format(
        script, branch[:], repos[:], current_time
    )
    cube.attributes[
        "Input"
    ] = "GEBCO_v{}_NO_REPEAT_LAT_LON.nc, created by  MAKE_GEBCO_CUBE.py, and the AMM15 unrotated coordinates file".format(YEAR)
    cube.attributes["Python version"] = platform.python_version()
    cube.attributes["System"] = platform.system()
    cube.attributes["Release"] = platform.release()


#%%
parser = argparse.ArgumentParser(description="Process inputs and output file paths")
parser.add_argument(
    "-a",
    "--AMM15_PATH",
    metavar="AMM15_cube_file",
    nargs=1,
    help="File location of AMM15 rotated grid",
    required=True,
)
parser.add_argument(
    '-i',
    '--INLSM_DIR',
    metavar='INLSM_DIR',
    nargs=1,
    help='Path to source LSM files', 
    required=True,
)

parser.add_argument(
    "-c",
    "--INCUBE_DIR",
    metavar="INCUBE_DIR",
    nargs=1,
    help="Path to source cube directory ",
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
    "-y",
    "--YEAR",
    metavar="YEAR",
    nargs=1,
    help="YEAR version of GEBCO e.g. 2020,2021 ",
    required=True,
)

args = parser.parse_args()
if not all([args.AMM15_PATH,args.INLSM_DIR, args.INCUBE_DIR, args.OUT_DIR,args.YEAR]):
    print(" Sorry All Arguments are required")
    sys.exit("Sorry all the arguments are required")


print(args.AMM15_PATH[0])
print("\n----------------------------------------------------\n")
print("Thanks, you have chosen: \n ")
print("      AMM15 grid file as {}\n".format(args.AMM15_PATH[0]))

if Path(args.AMM15_PATH[0]).is_file():
    print(" and the file {} exists.".format(args.AMM15_PATH[0]))
else:
    sys.exit("However, {} does not exist  so we exit here".format(args.AMM15_PATH[0]))
print("\n     the INPUT directory for LSM as {}\n".format(  args.INLSM_DIR[0] ))
if (Path(args.INLSM_DIR[0])).is_dir():
       print(" and the {} exists.".format (args.INLSM_DIR[0]))
else:
      sys.exit("However, {} does not exist  so we exit here".format(args.INLSM_DIR[0]))
print("\n     the INPUT directory for cube as {}\n".format(args.INCUBE_DIR[0]))
if (Path(args.INCUBE_DIR[0])).is_dir():
    print(" and the {} exists.".format(args.INCUBE_DIR[0]))
else:
    sys.exit("However, {} does not exist  so we exit here".format(args.INSUBE_DIR[0]))
print("\n     the OUTPUT directory for processed bathy as {}\n".format(args.OUT_DIR[0]))
if (Path(args.OUT_DIR[0])).is_dir():
    print(" and the {} exists.".format(args.OUT_DIR[0]))
else:
    sys.exit("However, {} does not exist  so we exit here".format(args.OUT_DIR[0]))
print("\n     and the GEBCO VERSION:  {}\n".format(args.YEAR[0]))


print("\n----------------------------------------------------\n")


print("\n Loading the AMM15  Cube to get base grid:\n")

AMM15_cube = iris.load(args.AMM15_PATH)[0]
X_LAT = np.array((AMM15_cube.coord("grid_latitude").points[:]))
X_LON = np.array((AMM15_cube.coord("grid_longitude").points[:]))


DEL_LAT = X_LAT[101] - X_LAT[100]
DEL_LON = X_LON[101] - X_LON[100]

EXTRA_LAT = int(100)
EXTRA_LON = int(100)

# Start,stop, increment
START_LAT = X_LAT[0] - EXTRA_LAT * DEL_LAT
START_LON = X_LON[0] - EXTRA_LON * DEL_LON

END_LAT = START_LAT + DEL_LAT * (np.size(X_LAT[:]) + 2 * EXTRA_LAT)
END_LON = START_LON + DEL_LON * (np.size(X_LON[:]) + 2 * EXTRA_LON)

inflate_lat = np.arange(START_LAT, END_LAT, DEL_LAT)
inflate_lon = np.arange(START_LON, END_LON, DEL_LON)

n_lat = np.size(inflate_lat)
n_lon = np.size(inflate_lon)


# Make Expansions

rotated_cs = iris.coord_systems.RotatedGeogCS(37.5, 177.5)

# set up rotated cube using AMM15 read in lats and lons

grid_latitude = DimCoord(
    inflate_lat, standard_name="grid_latitude", units="degrees", coord_system=rotated_cs
)

grid_longitude = DimCoord(
    inflate_lon,
    standard_name="grid_longitude",
    units="degrees",
    coord_system=rotated_cs,
)

print("\n Creating the new expanded cube:\n")

expandcube = Cube(
    np.zeros((n_lat, n_lon), np.float32),
    dim_coords_and_dims=[(grid_latitude, 0), (grid_longitude, 1)],
    standard_name="sea_floor_depth_below_geoid",
)


print("\n Saving the expanded domain to expand_test:\n")

iris.save(expandcube, "{}/expand_test.nc".format(args.OUT_DIR[0]))


# we want to save the 2d rotated lats and lons for the extended grid in NEMO format


output_nemo_coords(inflate_lat, inflate_lon, args.OUT_DIR[0])


# Read in the cube of gebco data this was created by MAKE_GEBCO_CUBE.py

#GEBCO_RAW_cube = iris.load(
#    "{}/GEBCO_v2020_NO_REPEAT_LAT_LON.nc".format(args.INCUBE_DIR[0])
#)[0]

GEBCO_RAW_cube = iris.load('{}/GEBCO_CUBE.nc'.format( args.INCUBE_DIR[0] ))[0]


## Read in the in the LSM for the GEBCO data
YEAR=args.YEAR[0]
GEBCO_LSM_fp = Dataset('{}/NWS_CUT_GEBCO_{}_TID.nc'.format( args.INLSM_DIR[0],YEAR ),'r')
GEBCO_LSM = GEBCO_LSM_fp.variables['tid'][:]


#GEBCO_RAW_cube.data[np.where(GEBCO_LSM.data[:] ==0) ] = np.nan


expandcube.coord_system = AMM15_cube.coord_system

# GEBCO_RAW_cube.data = GEBCO_RAW_cube.lazy_data().rechunk([10,None])
print(GEBCO_RAW_cube.lazy_data().chunks)


print(" Size of domain we want to split is : ", np.size(expandcube.data[:, 0]))

# Sub size is a small patch to interpolate to to keep the memory usage down of iris regrid
# we make make the src grid bound this plus some extra points for consistent infilling
# make smaller for smaller RAM
# SUB_SIZE  = 100
SUB_SIZE = 100

div, remain = np.divmod(np.size(expandcube.data[:, 0]), SUB_SIZE)
print(
    " This means we need {}  domains of size {} and a remainder domain of size {}".format(
        div, SUB_SIZE, remain
    )
)
SECTION = 0

os.system("mkdir  -p {}/SUBSECTION/".format(args.OUT_DIR[0]))
# Clear out incase there is anything there
os.system("rm {}/SUBSECTION/*".format(args.OUT_DIR[0]))
#%%
while SECTION < div + 1:

    print("Doing SECTION {} of {} ".format(SECTION, div))
    if SECTION < div:
        subsection_cube = expandcube[SECTION * SUB_SIZE : (SECTION + 1) * SUB_SIZE, :]
    if SECTION == div:
        subsection_cube = expandcube[SECTION * SUB_SIZE :, :]

    X, Y = np.meshgrid(
        subsection_cube.coord("grid_longitude").points,
        subsection_cube.coord("grid_latitude").points,
    )

    lons, lats = iris.analysis.cartography.unrotate_pole(X, Y, 177.5, 37.5)

    # find closest indices
    # Min
    values, lat_idx = find_nearest(
        GEBCO_RAW_cube.coord("latitude").points[:], np.min(lats)
    )
    lat_min_idx = lat_idx - 500  # ( for safety for fill interp)

    values, lon_idx = find_nearest(
        GEBCO_RAW_cube.coord("longitude").points[:], np.min(lons)
    )
    lon_min_idx = lon_idx - 500  # ( for safety for fill interp)

    # Max
    values, lat_idx = find_nearest(
        GEBCO_RAW_cube.coord("latitude").points[:], np.max(lats)
    )
    lat_max_idx = lat_idx + 500  # (for safety for fill interp)

    values, lon_idx = find_nearest(
        GEBCO_RAW_cube.coord("longitude").points[:], np.max(lons)
    )
    lon_max_idx = lon_idx + 500  # ( for safety for fill interp)

    print("\n SECTION {} of {}\n".format(SECTION, div))

    subsection_gebco_raw_cube = GEBCO_RAW_cube[
        lat_min_idx:lat_max_idx, lon_min_idx:lon_max_idx]

    subsection_gebco_lsm = GEBCO_LSM[
        lat_min_idx:lat_max_idx, lon_min_idx:lon_max_idx]

    subsection_gebco_raw_cube.data[np.where(subsection_gebco_lsm.data[:] ==0) ] = np.nan

    print("About to save and regrid")
    set_history(subsection_gebco_raw_cube)
    iris.save(
        subsection_gebco_raw_cube.regrid(
            subsection_cube, iris.analysis.Linear(extrapolation_mode="mask")
        ),
        "{}/SUBSECTION/{:05d}_SUBSECTION_CUBE.nc".format(args.OUT_DIR[0], SECTION),
    )

    print("About to Nan")
    with ProgressBar():
        subsection_gebco_raw_cube.data[
            subsection_gebco_raw_cube.data > 1.0e33
        ] = np.nan

    print("About to Fill")
    with ProgressBar():
        subsection_gebco_raw_cube.data = fill(subsection_gebco_raw_cube.data)

    print("About to save and regrid")
    set_history(subsection_gebco_raw_cube)
    format(args.OUT_DIR[0])
    with ProgressBar():
        iris.save(
            subsection_gebco_raw_cube.regrid(
                subsection_cube, iris.analysis.Linear(extrapolation_mode="extrapolate")
            ),
            "{}/SUBSECTION/{:05d}_GEBCO_FILL_SUBSECTION_CUBE.nc".format(
                args.OUT_DIR[0], SECTION
            ),
        )
    print("SECTION {} of {} done".format(SECTION, div))

    SECTION += 1

print("all sections done")
#%%
print("{}/SUBSECTION/00???_GEBCO_FILL_SUBSECTION_CUBE.nc".format(args.OUT_DIR[0]))
print("Files to stitch are:")
for file_name in glob.iglob(
    "{}/SUBSECTION/00???_GEBCO_FILL_SUBSECTION_CUBE.nc".format(args.OUT_DIR[0]),
    recursive=True,
):
    print(file_name)
#%%
EXTRAPOLATE_GEBCO_ON_EXPANDAMM15 = xr.open_mfdataset(
    "{}/SUBSECTION/00???_GEBCO_FILL_SUBSECTION_CUBE.nc".format(args.OUT_DIR[0]),
    combine="nested",
    concat_dim="grid_latitude",
    parallel=True,
)

MASK_GEBCO_ON_EXPANDAMM15 = xr.open_mfdataset(
    "{}/SUBSECTION/00???_SUBSECTION_CUBE.nc".format(args.OUT_DIR[0]),
    combine="nested",
    concat_dim="grid_latitude",
    parallel=True,
)

EXTRA_LAT = int(100)
EXTRA_LON = int(100)

MASK_EXTRAPOLATE = MASK_GEBCO_ON_EXPANDAMM15.copy()
MASK_EXTRAPOLATE.sea_floor_depth_below_geoid[
    EXTRA_LAT:-EXTRA_LAT, EXTRA_LON:-EXTRA_LON
] = EXTRAPOLATE_GEBCO_ON_EXPANDAMM15.sea_floor_depth_below_geoid[
    EXTRA_LAT:-EXTRA_LAT, EXTRA_LON:-EXTRA_LON
]


with ProgressBar():
    MASK_GEBCO_ON_EXPANDAMM15.to_netcdf(
        "{}/MASK_GEBCO_vDec{}_ON_EXPAND_AMM15.nc".format(args.OUT_DIR[0],YEAR)
    )
with ProgressBar():
    EXTRAPOLATE_GEBCO_ON_EXPANDAMM15.to_netcdf(
        "{}/EXTRAPOLATE_GEBCO_vDec{}_ON_EXPAND_AMM15.nc".format(args.OUT_DIR[0],YEAR)
    )
with ProgressBar():
    MASK_EXTRAPOLATE.to_netcdf(
        "{}/MASK_EXTRAPOLATE_GEBCO_vDec{}_ON_EXPAND_AMM15.nc".format(
            args.OUT_DIR[0],YEAR
        )
    )
