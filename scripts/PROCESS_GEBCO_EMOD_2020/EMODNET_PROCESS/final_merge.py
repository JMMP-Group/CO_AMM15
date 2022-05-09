import xarray as xr
import numpy as np

import xarray as xr


import sys
PATH = sys.argv[1]


dsa = xr.open_mfdataset('%s/[1234]merge.nc'%(PATH),
                         combine = 'nested',  
                         concat_dim = 'lat', 
                         parallel=True,
                         chunks=({"lat": 3000, "lon": -1}))


from dask.diagnostics import ProgressBar

with ProgressBar():
  dsa.elevation.to_netcdf('%s/ALLmerge.nc'%(PATH))


print('All Done')

