import xarray as xr
import numpy as np

import xarray as xr


import sys
FILE = sys.argv[1]
NO = sys.argv[2]

print('%s*.nc'%(FILE))

## as there is overlap we combine by nested and specify long as concat dimension later we remove overlap
dsa = xr.open_mfdataset('%s*.nc'%(FILE), combine = 'nested', concat_dim = 'lon', parallel=True)


from dask.diagnostics import ProgressBar

print('%smerge.nc'%(NO))
with ProgressBar():
  dsa.elevation.to_netcdf('%smerge.nc'%(NO))


print('All Done')

