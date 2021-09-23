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
HPGEfile = '/data/users/dbruciaf/JMMP_HPGE/amm15/MEs_2env_0.1_0.07/maximum_hpge_ref.nc'
RMAXfile = '/data/users/dbruciaf/JMMP_HPGE/amm15/tmp_jelt_diego/bathymetry.amm15_MEs_2env_0.1_0.07.nc'

trshl = 0.08

ds_pe = xr.open_dataset(HPGEfile)
hpge  = ds_pe['max_hpge_1']
ds_rr = xr.open_dataset(RMAXfile)
rmaxr = ds_rr['rmax0_2']
env_r = ds_rr['hbatt_2']

fig, ax = plt.subplots(ncols=1, nrows=1)
ax.scatter(hpge.data.flatten(),rmaxr.data.flatten(),s=20)
ax.plot([trshl,trshl],[0.,rmaxr.max()])
ax.set_xlabel('max hpge [$m\;s^{-1}$]')
ax.set_ylabel('rmax')
plt.show()
