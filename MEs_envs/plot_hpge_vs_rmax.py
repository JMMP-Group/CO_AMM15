#!/usr/bin/env python

#     |------------------------------------------------------------|
#     | Author: Diego Bruciaferri                                  |
#     | Date and place: 07-09-2021, Met Office, UK                 |
#     |------------------------------------------------------------|


import os
from os.path import join, isfile, basename, splitext
import glob
import numpy as np
from matplotlib import pyplot as plt
import xarray as xr

# ==============================================================================
# Input files
# ==============================================================================

# Path of the file containing max HPGE spurious currents velocity field
HPGEfile = '/data/users/dbruciaf/mod_dev/pgm18/vgrid/hpge_MEs_2env_0.3_0.1/maximum_hpge.nc'
RMAXfile = '/data/users/dbruciaf/mod_dev/files4bdy_gulf18/bathymetry.gulf18_MEs_2env_0.3_0.1.nc'

#trshl = 0.06

ds_pe = xr.open_dataset(HPGEfile)
hpge  = ds_pe['max_hpge_0']
ds_rr = xr.open_dataset(RMAXfile)
rmaxr = ds_rr['rmax0_1']
env_r = ds_rr['hbatt_1']
ds_ro = xr.open_dataset(OPTIfile)
env_o = ds_ro['hbatt_1']
bathy = ds_ro['Bathymetry']

fig, ax = plt.subplots(ncols=1, nrows=1)
ax.scatter(hpge.data.flatten(),rmaxr.data.flatten(),s=20)
ax.plot([trshl,trshl],[0.,rmaxr.max()])
ax.set_xlabel('max hpge [$m\;s^{-1}$]')
ax.set_ylabel('rmax')
plt.show()
