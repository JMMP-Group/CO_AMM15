#!/usr/bin/env python

#     |------------------------------------------------------------|
#     | This module creates general envelope surfaces to be        |
#     | used to generate a Multi-Envelope s-coordinate             |
#     | vertical grid.                                             |
#     |                                                            |
#     |                                                            |
#     |                      -o#&&*''''?d:>b\_                     |
#     |                 _o/"`''  '',, dMF9MMMMMHo_                 |
#     |              .o&#'        `"MbHMMMMMMMMMMMHo.              |
#     |            .o"" '         vodM*$&&HMMMMMMMMMM?.            |
#     |           ,'              $M&ood,~'`(&##MMMMMMH\           |
#     |          /               ,MMMMMMM#b?#bobMMMMHMMML          |
#     |         &              ?MMMMMMMMMMMMMMMMM7MMM$R*Hk         |
#     |        ?$.            :MMMMMMMMMMMMMMMMMMM/HMMM|`*L        |
#     |       |               |MMMMMMMMMMMMMMMMMMMMbMH'   T,       |
#     |       $H#:            `*MMMMMMMMMMMMMMMMMMMMb#}'  `?       |
#     |       ]MMH#             ""*""""*#MMMMMMMMMMMMM'    -       |
#     |       MMMMMb_                   |MMMMMMMMMMMP'     :       |
#     |       HMMMMMMMHo                 `MMMMMMMMMT       .       |
#     |       ?MMMMMMMMP                  9MMMMMMMM}       -       |
#     |       -?MMMMMMM                  |MMMMMMMMM?,d-    '       |
#     |        :|MMMMMM-                 `MMMMMMMT .M|.   :        |
#     |         .9MMM[                    &MMMMM*' `'    .         |
#     |          :9MMk                    `MMM#"        -          |
#     |            &M}                     `          .-           |
#     |             `&.                             .              |
#     |               `~,   .                     ./               |
#     |                   . _                  .-                  |
#     |                     '`--._,dd###pp=""'                     |
#     |                                                            |
#     |                                                            |
#     | Author: Diego Bruciaferri                                  |
#     | Date and place: 26-07-2021, Met Office, UK                 |
#     |------------------------------------------------------------|


import sys
import os
from os.path import isfile, basename, splitext
import numpy as np
from matplotlib import pyplot as plt
import xarray as xr

from lib_env import *

# ==============================================================================
# 1. Checking for input files
# ==============================================================================

msg_info('GRID ENVELOPES GENERATOR', main=True)

# Checking input file for envelopes
if len(sys.argv) != 2:
   raise TypeError(
         'You need to specify the absolute path of the envelopes input file'
   )

# Reading envelopes infos
env_file = sys.argv[1]
msg_info('Reading envelopes parameters ....')
envInfo = read_envInfo(env_file)

# Reading bathymetry and horizontal grid
msg_info('Reading bathymetry data ... ')
ds_bathy = xr.open_dataset(envInfo.bathyFile).squeeze()
bathy = ds_bathy["Bathymetry"].squeeze()

msg_info('Reading horizontal grid data ... ')
ds_grid = xr.open_dataset(envInfo.hgridFile) 
glamt = ds_grid["glamt"]
gphit = ds_grid["gphit"]

ds_env = ds_bathy.copy()

# Defining local variables -----------------------------------------------------

num_env = len(envInfo.e_min_ofs)
ni      = glamt.shape[1]
nj      = glamt.shape[0]

e_min_ofs = envInfo.e_min_ofs
e_max_dep = envInfo.e_max_dep
e_loc_vel = envInfo.e_loc_vel
e_loc_var = envInfo.e_loc_var
e_loc_vmx = envInfo.e_loc_vmx
e_loc_rmx = envInfo.e_loc_rmx
e_loc_hal = envInfo.e_loc_hal
e_glo_rmx = envInfo.e_glo_rmx

# Computing LSM
lsm = xr.where(bathy > 0, 1, 0) 

#--------------------------------------------------------------------------------------
# Cretaing envelopes
for env in range(num_env):

    msg = 'ENVELOPE ' + str(env) + ': min_ofs = '+ str(e_min_ofs[env]) + \
          ', max_dep = ' + str(e_max_dep[env]) + ', glo_rmx = ' + str(e_glo_rmx[env])
    msg_info(msg, main=True)

    env_bathy = bathy.copy()

    # =============================================================
    # 1. GEOMETRY of the envelope
    #
    # e_min_ofs[env] > 0: offset to be added to surf to comupute 
    #                     minimum depth of envelope env
    # e_min_ofs[env] < 0: e_min_ofs[env] is equal to e_max_dep[env]
    # e_max_dep[env] < 0: e_max_dep[env] = np.amax(bathy)

    msg = '1. Computing the geometry of the envelope'
    msg_info(msg)
    
    if env == 0:
       # For the first envelope, the upper envelope 
       # is the ocean surface at rest
       surf = xr.full_like(bathy,0.,dtype=np.double)
    else:
       # For deeper envelopes, the upper 
       # envelope is the previous envelope
       surf = ds_env["hbatt_"+str(env)]

    if e_max_dep[env] < 0.: e_max_dep[env] = np.nanmax(bathy)

    # CASE A flat envelope
    if e_min_ofs[env] < 0.:                     
       msg = 'Generating a flat envelope'
       msg_info(msg)
       hbatt = xr.full_like(bathy,e_max_dep[env],dtype=np.double)
       
    # CASE B general envelope
    else:
       hbatt = calc_zenv(env_bathy, surf, e_min_ofs[env], e_max_dep[env])  
 
    # For the first envelope, we follow NEMO approach
    #if env == 0:
    #   hbatt *= lsm

    hbatt.plot.pcolormesh(add_colorbar=True, add_labels=True, \
                          cbar_kwargs=dict(pad=0.15, shrink=1, \
                          label='Envelope '+ str(env)))

    plt.show()
    # =============================================================
    # 2. SMOOTHING the envelope
  
    msg = '2. Smoothing the envelope'
    msg_info(msg)  

    # Computing then MB06 Slope Parameter for the raw envelope
    rmax_raw = np.amax(calc_rmax(hbatt)*lsm).values
    msg = '   Slope parameter of raw envelope: rmax = ' + str(rmax_raw)
    msg_info(msg)

    if len(e_loc_vel[env]) == 0:
       hbatt_smt = hbatt.copy()
    else:
       # LOCAL SMOOTHING
       hal = e_loc_hal[env]
       r0x = e_loc_rmx[env]
       msg = ("   Local smoothing of areas of the raw" +\
              " envelope where HPGE > " + str(e_loc_vmx[env]) + " m/s\n"+\
              "   Envelope  : " + str(env) + "\n"+\
              "   Local rmax: " + str(r0x) + "\n"+\
              "   Local halo: " + str(hal) + " cells")
       msg_info(msg)
              
       env_tmp = hbatt.copy()
       msk_pge = env_tmp.copy()*0.
 
       for m in range(len(e_loc_vel[env])):

           filename = e_loc_vel[env][m]
           varname  = e_loc_var[env][m]
           ds_vel = xr.open_dataset(filename)
           hpge = ds_vel[varname].data

           nj = hpge.shape[0]
           ni = hpge.shape[1]
           for j in range(nj):
               for i in range(ni):
                   max_hpge = hpge[j,i]
                   if max_hpge >= e_loc_vmx[env]:
                      msk_pge[j-hal:j+hal+1,i-hal:i+hal+1] = 1
                          
       msg = '   Total number of points with HPGE > ' + str(e_loc_vmx[env]) + ' m/s: ' + str(np.nansum(msk_pge)) 
       msg_info(msg,)

       msk_pge.plot.pcolormesh(add_colorbar=True, add_labels=True, \
                               cbar_kwargs=dict(pad=0.15, shrink=1, \
                               label='HPGE smoothing mask'))
       plt.show()

       # smoothing with Martinho & Batteen 2006 
       hbatt_smt = smooth_MB06(env_tmp, e_loc_rmx[env])
       
       # applying smoothing only where HPGE are large
       WRK = hbatt_smt.data
       TMP = env_tmp.data
       WRK[msk_pge==0] = TMP[msk_pge==0]
       hbatt_smt.data = WRK

    # GLOBAL SMOOTHING
    if e_glo_rmx[env] > 0:

       msg = ('   Global smoothing of the raw envelope')
       msg_info(msg)

       #da_wrk = hbatt_smt.copy()
       #env_wrk = np.copy(hbatt_smt.data)

       if env == 0:
          # for the first envelope, we follow NEMO approach
          # set first land point adjacent to a wet cell to
          # min_dep as this needs to be included in smoothing
          cst_lsm = lsm.rolling({dim: 3 for dim in lsm.dims}, min_periods=2).sum()
          cst_lsm = cst_lsm.shift({dim: -1 for dim in lsm.dims})
          cst_lsm = (cst_lsm > 0) & (lsm == 0)
          da_wrk = hbatt_smt.where(cst_lsm == 0, e_min_ofs[env])
       else:
          da_wrk = hbatt_smt.copy()

       # smoothing with Martinho & Batteen 2006 
       hbatt_smt = smooth_MB06(da_wrk, e_glo_rmx[env])#, tol=2.5e-8)

    # Computing then MB06 Slope Parameter for the smoothed envelope
    rmax0_smt = calc_rmax(hbatt_smt)*lsm
    rmax0_smt.plot.pcolormesh(add_colorbar=True, add_labels=True, \
                              cbar_kwargs=dict(pad=0.15, shrink=1, \
                              label='Slope parameter'))
    plt.show()
    rmax_smt = np.amax(rmax0_smt).data
    msg = '   Slope parameter of smoothed envelope: rmax = ' + str(rmax_smt)
    msg_info(msg)

    # Saving envelope DataArray
    ds_env["hbatt_"+str(env+1)] = hbatt_smt
    ds_env["rmax0_"+str(env+1)] = rmax0_smt

# -------------------------------------------------------------------------------------   
# Writing the bathy_meter.nc file

msg = 'WRITING the bathy_meter.nc FILE'
msg_info(msg)

out_name = splitext(basename(env_file))[0] 

out_file = "bathymetry." + out_name + ".nc"

ds_env.to_netcdf(out_file)

