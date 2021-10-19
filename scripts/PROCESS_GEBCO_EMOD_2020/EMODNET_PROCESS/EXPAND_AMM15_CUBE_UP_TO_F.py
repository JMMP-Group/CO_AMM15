""" Interpolate EMODNET data onto Expanded AMM15 grid
We read in the rotated amm15 grid. We create a grid based on this but that is wider
so that later smoothing will not effect the bdy regions.
Here is it set arbitarliy as a 100 points extra NSWE(rotated space)

If we were just considering the AMM15 we would enforce the existing operational LSM
on this domain in the next step as well as an LAT correction derived from Surge models

However, as we go beyond the AMM15 edges we need to use a LSM as derived from EMODNET itself
We can extrapolate bathy values inside the AMM15 so that we get values in W&D areas,
when we later impose the AMM15 mask.

Thus we have a version of the data that is extrapolated into the land
and a version that has the EMODET LSM applied directly.

Later we use the EMODET LSM in the areas beyond the AMM15 extents and the oper LSM in the 
inner true AMM15 region

Note on Spice we can use enough RAM to process the iris regird all in one but on
more reasonable machsin etis is not possible and iris regirds cannot chunk in the horizontal

So we have a rok aroudn where we maek a small section of destination grid
and a small section of src data that covers the destination gris plus extra (1000)
so the infilling is done consistently

then at the end we can rejoin all the interpolated subsections into a single file


Parameters
----------
AMM15_PATH : str
    The location of the rotated AMM15 grid
INCUBE_DIR : str
    The path the input cube of bathymetry
OUT_DIR : str
    The path write the output cube to

Returns
-------
EMODNET interpolated on the expanded AMM15 grid

"""
#%%
import iris
from iris.cube import Cube
from iris.coords import DimCoord
from dask.diagnostics import ProgressBar

import numpy as np

import os, glob


import argparse
from pathlib import Path

import sys, platform

sys.path.append(r'../')
from mask_tools import fill

import xarray as xr

from datetime import datetime
import subprocess



def find_nearest(array, value):
   ''' 
   Finds the nearest value and index of that value in an ordered array to a specifed value

            Parameters:
                    array[] (float): ordered array of values (lats/lons)
                    value (float): Desired point 

            Returns:
                    array[idx] (float): value at array nearest to value
                    idx(int) : index of array closes to value
                    
   '''

   array = np.asarray(array)
   idx = (np.abs(array - value)).argmin()
   return array[idx],idx

def set_history(cube):
    '''
    Adds a record to the cubes history of when and how it was made 
   
            Parameters:
                    cube(iris cube) : cube to be stored in netcdf
    '''
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
    cube.attributes[ 'History' ] = "Created by {} from branch {} of {} on {} ".format(script,branch[:],repos[:],current_time)
    cube.attributes[ 'Input' ]   = "EMODNET_CUBE.nc, created by  MAKE_EMODNET_CUBE.py,  and the AMM15 unrotated coordinates file"
    cube.attributes[ 'Python version' ] = platform.python_version()
    cube.attributes[ 'System' ]  = platform.system()
    cube.attributes[ 'Release' ] = platform.release()




#%%
parser = argparse.ArgumentParser(description='Process inputs and output file paths')
parser.add_argument('-a','--AMM15_PATH', metavar='AMM15_cube_file', nargs=1,
                    help='File location of AMM15 rotated grid', required=True)
parser.add_argument('-c','--INCUBE_DIR',metavar='INCUBE_DIR',  nargs=1,
                    help='Path to source cube directory ', required=True)
parser.add_argument( '-o','--OUT_DIR',metavar='OUT_DIR', nargs=1,
                    help='Path to output to ', required=True) 

args = parser.parse_args()
if not all([args.AMM15_PATH,  args.INCUBE_DIR, args.OUT_DIR]):
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

print("\n----------------------------------------------------\n")


print("\n Loading the AMM15  Cube to get base grid:\n")

AMM15_cube = iris.load(args.AMM15_PATH)[0]
X_LAT = np.array((AMM15_cube.coord('grid_latitude').points[:]))
X_LON = np.array((AMM15_cube.coord('grid_longitude').points[:]))


DEL_LAT=X_LAT[101]-X_LAT[100]
DEL_LON=X_LON[101]-X_LON[100]

extra_lat=int(100)
extra_lon=int(100)

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


# Read in the cube of gebco data this was created by MAKE_EMODNET_CUBE.py

EMODNET_RAW_cube = iris.load('{}/EMODNET_v2020_NO_REPEAT_LAT_LON.nc'.format( args.INCUBE_DIR[0] ))[0]


expandcube.coord_system = AMM15_cube.coord_system

# EMODNET_RAW_cube.data = EMODNET_RAW_cube.lazy_data().rechunk([10,None])
print(EMODNET_RAW_cube.lazy_data().chunks)


print ( " Size of domain we want to split is : ", np.size(expandcube.data[:,0]) )

# Sub size is a small patch to interpolate to to keep the memory usage down of iris regrid
# we make make the src grid bound this plus some extra points for consistent infilling
# make smaller for smaller RAM
#sub_size  = 100
sub_size  = 300

div, remain = np.divmod( np.size(expandcube.data[ : , 0]), sub_size )
print ( " This means we need {}  domains of size {} and a remainder domain of size {}".format(div,sub_size,remain) )
section = 0

os.system("mkdir  -p {}/SUBSECTION/".format(args.OUT_DIR[0]) )
# Clear out incase there is anything there
os.system("rm {}/SUBSECTION/*".format(args.OUT_DIR[0]) )
#%%
while section  < div+1:
  
  print( "Doing Section {} of {} ".format(section,div))
  if(section < div):
     subsection_cube = expandcube[ section*sub_size:(section+1)*sub_size,:]
  if(section == div):
     subsection_cube = expandcube[ section*sub_size:,:]

  X, Y = np.meshgrid(subsection_cube.coord('grid_longitude').points,subsection_cube.coord('grid_latitude').points)

  lons, lats = iris.analysis.cartography.unrotate_pole(X ,Y    , 177.5, 37.5)

#find closest indices
# Min
  value,lat_idx = find_nearest(EMODNET_RAW_cube.coord('latitude').points[:], np.min(lats ))
  lat_min_idx = lat_idx - 500 # ( for safety for fill interp)
   

  value,lon_idx = find_nearest(EMODNET_RAW_cube.coord('longitude').points[:], np.min(lons ))
  lon_min_idx = lon_idx - 500 # ( for safety for fill interp)
    


# Max
  value,lat_idx = find_nearest(EMODNET_RAW_cube.coord('latitude').points[:], np.max(lats ))
  lat_max_idx = lat_idx + 500 # (for safety for fill interp)
    

  value,lon_idx = find_nearest(EMODNET_RAW_cube.coord('longitude').points[:], np.max(lons ))
  lon_max_idx = lon_idx + 500 # ( for safety for fill interp)
    



  print("\n section {} of {}\n".format(section,div) )
  subsection_emodnet_raw_cube = EMODNET_RAW_cube[lat_min_idx:lat_max_idx, lon_min_idx:lon_max_idx ]




  print("About to save and regrid")
  set_history(subsection_emodnet_raw_cube )
  iris.save( subsection_emodnet_raw_cube.regrid(subsection_cube, iris.analysis.Linear(extrapolation_mode='mask') ), '{}/SUBSECTION/{:05d}_SUBSECTION_CUBE.nc'.format(args.OUT_DIR[0],section))


  print("About to Nan")
  with ProgressBar():
     subsection_emodnet_raw_cube.data[subsection_emodnet_raw_cube.data > 1.e33 ] = np.nan

  print("About to Fill")
  with ProgressBar():
     subsection_emodnet_raw_cube.data = fill(subsection_emodnet_raw_cube.data)

  print("About to save and regrid")
  set_history(subsection_emodnet_raw_cube )
  format(  args.OUT_DIR[0] )
  with ProgressBar():
     iris.save(subsection_emodnet_raw_cube.regrid(subsection_cube, iris.analysis.Linear(extrapolation_mode='extrapolate')) ,'{}/SUBSECTION/{:05d}_FILL_SUBSECTION_CUBE.nc'.format(args.OUT_DIR[0],section))
  print("Section {} of {} done".format(section,div))


  section += 1

print("all sections done")
#%%
print('{}/SUBSECTION/00???_FILL_SUBSECTION_CUBE.nc'.format(args.OUT_DIR[0]))
print("Files to stitch are:")
for file_name in glob.iglob('{}/SUBSECTION/00???_FILL_SUBSECTION_CUBE.nc'.format(args.OUT_DIR[0]), recursive=True):
  print(file_name)
#%%
EXTRAPOLATE_EMODNET_ON_EXPANDAMM15 = xr.open_mfdataset(
        '{}/SUBSECTION/00???_FILL_SUBSECTION_CUBE.nc'.format(args.OUT_DIR[0]), 
        combine = 'nested', 
        concat_dim = 'grid_latitude', 
        parallel=True)

MASK_EMODNET_ON_EXPANDAMM15 = xr.open_mfdataset(
        '{}/SUBSECTION/00???_SUBSECTION_CUBE.nc'.format(args.OUT_DIR[0]),
        combine = 'nested', 
        concat_dim = 'grid_latitude', 
        parallel=True)

extra_lat=int(100)
extra_lon=int(100)

MASK_EXTRAPOLATE = MASK_EMODNET_ON_EXPANDAMM15.copy()
MASK_EXTRAPOLATE.sea_floor_depth_below_geoid[extra_lat:-extra_lat,extra_lon:-extra_lon] = (
EXTRAPOLATE_EMODNET_ON_EXPANDAMM15.sea_floor_depth_below_geoid [extra_lat:-extra_lat,extra_lon:-extra_lon]   )



with ProgressBar():
  MASK_EMODNET_ON_EXPANDAMM15.to_netcdf('{}/MASK_EMODNET_vDec2020_ON_EXPAND_AMM15.nc'.format(args.OUT_DIR[0]))
with ProgressBar():
  EXTRAPOLATE_EMODNET_ON_EXPANDAMM15.to_netcdf('{}/EXTRAPOLATE_EMODNET_vDec2020_ON_EXPAND_AMM15.nc'.format(args.OUT_DIR[0]))
with ProgressBar():
  MASK_EXTRAPOLATE.to_netcdf('{}/MASK_EXTRAPOLATE_EMODNET_vDec2020_ON_EXPAND_AMM15.nc'.format(args.OUT_DIR[0]))

