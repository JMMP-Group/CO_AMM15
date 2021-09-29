import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

import xarray as xr


import sys
FILE = sys.argv[1]
NO = sys.argv[2]

dsa = xr.open_mfdataset('../NC/%s*.nc'%(FILE), parallel=True)


from dask.diagnostics import ProgressBar

with ProgressBar():
  dsa.elevation.to_netcdf('../MERGE/%smerge.nc'%(NO))


print('All Done')

