#!/usr/bin/env python

from typing import Tuple
import numpy as np
import netCDF4 as nc4
import xarray as xr
from xarray import Dataset, DataArray 

#=======================================================================================
def read_envInfo(filepath):
    '''
    This function reads the parameters which will be used
    to generate the model grid envelopes.

    If the file is not in the right format and does not contain
    all the needed infos an error message is returned.

    filepath: string
    '''

    import importlib.machinery as imp

    loader  = imp.SourceFileLoader(filepath,filepath)
    envInfo = loader.load_module()

    attr = ['bathyFile', 'hgridFile', \
            'e_min_ofs', 'e_max_dep', \
            'e_loc_vel', 'e_loc_var', 'e_loc_vmx', \
            'e_loc_rmx', 'e_loc_hal', 'e_glo_rmx']

    errmsg = False
    for a in attr:
        if not hasattr(envInfo,a):
           errmsg = True
           attrerr = a
        else:
           if getattr(envInfo,a) == "":
              errmsg = True
              attrerr = a

    if errmsg: 
       raise AttributeError(
             'Attribute ' + attrerr + ' is missing'
       )

    return envInfo

#=======================================================================================
def calc_zenv(bathy, surf, offset, max_dep):
    """
    Constructs an envelop bathymetry 
    for the sigma levels definition.
       
    INPUT:
    *) bathy: the bathymetry field which will be 
              used to compute the envelope.
    *) surf:  2D field used as upper surface to model 
              the envelope
    *) offset: offset to be added to surf to calc 
               the minimum depth of the envelope.

    *) max_depth : maximum depth of the envelope.
    """

    zenv   = bathy.copy()
    env_up = surf.copy()

    # Set maximum depth of the envelope
    zenv = np.minimum(zenv, max_dep)

    # Set minimum depth of the envelope
    env_up += offset
    zenv = np.maximum(zenv, env_up)

    return zenv

#=======================================================================================
def calc_rmax(depth: DataArray) -> DataArray:
    """
    Calculate rmax: measure of steepness
    This function returns the slope paramater field

    r = abs(Hb - Ha) / (Ha + Hb)

    where Ha and Hb are the depths of adjacent grid cells (Mellor et al 1998).

    Reference:
    *) Mellor, Oey & Ezer, J Atm. Oce. Tech. 15(5):1122-1131, 1998.

    Parameters
    ----------
    depth: DataArray
        Bottom depth (units: m).

    Returns
    -------
    DataArray
        2D slope parameter (units: None)

    Notes
    -----
    This function uses a "conservative approach" and rmax is overestimated.
    rmax at T points is the maximum rmax estimated at any adjacent U/V point.
    """
    # Mask land
    depth = depth.where(depth > 0)

    # Loop over x and y
    both_rmax = []
    for dim in depth.dims:

        # Compute rmax
        rolled = depth.rolling({dim: 2}).construct("window_dim")
        diff = rolled.diff("window_dim").squeeze("window_dim")
        rmax = np.abs(diff) / rolled.sum("window_dim")

        # Construct dimension with velocity points adjacent to any T point
        # We need to shift as we rolled twice
        rmax = rmax.rolling({dim: 2}).construct("vel_points")
        rmax = rmax.shift({dim: -1})

        both_rmax.append(rmax)

    # Find maximum rmax at adjacent U/V points
    rmax = xr.concat(both_rmax, "vel_points")
    rmax = rmax.max("vel_points", skipna=True)

    # Mask halo points
    for dim in rmax.dims:
        rmax[{dim: [0, -1]}] = 0

    return rmax.fillna(0)

#=======================================================================================
def smooth_MB06(
    depth: DataArray,
    rmax: float,
    tol: float = 1.0e-8,
    max_iter: int = 10_000,
) -> DataArray:
    """
    Direct iterative method of Martinho and Batteen (2006) consistent
    with NEMO implementation.

    The algorithm ensures that

                H_ij - H_n
                ---------- < rmax
                H_ij + H_n

    where H_ij is the depth at some point (i,j) and H_n is the
    neighbouring depth in the east, west, south or north direction.

    Reference:
    *) Martinho & Batteen, Oce. Mod. 13(2):166-175, 2006.

    Parameters
    ----------
    depth: DataArray
        Bottom depth.
    rmax: float
        Maximum slope parameter allowed
    tol: float, default = 1.0e-8
        Tolerance for the iterative method
    max_iter: int, default = 10000
        Maximum number of iterations

    Returns
    -------
    DataArray
        Smooth version of the bottom topography with
        a maximum slope parameter < rmax.
    """

    # Set scaling factor used for smoothing
    zrfact = (1.0 - rmax) / (1.0 + rmax)

    # Initialize envelope bathymetry
    zenv = depth

    for _ in range(max_iter):

        # Initialize lists of DataArrays to concatenate
        all_ztmp = []
        all_zr = []
        for dim in zenv.dims:

            # Shifted arrays
            zenv_m1 = zenv.shift({dim: -1})
            zenv_p1 = zenv.shift({dim: +1})

            # Compute zr
            zr = (zenv_m1 - zenv) / (zenv_m1 + zenv)
            zr = zr.where((zenv > 0) & (zenv_m1 > 0), 0)
            for dim_name in zenv.dims:
                zr[{dim_name: -1}] = 0
            all_zr += [zr]

            # Compute ztmp
            zr_p1 = zr.shift({dim: +1})
            all_ztmp += [zenv.where(zr <= rmax, zenv_m1 * zrfact)]
            all_ztmp += [zenv.where(zr_p1 >= -rmax, zenv_p1 * zrfact)]

        # Update envelope bathymetry
        zenv = xr.concat([zenv] + all_ztmp, "dummy_dim").max("dummy_dim")

        # Check target rmax
        zr = xr.concat(all_zr, "dummy_dim")
        #print(np.nanmax(np.abs(zr)))
        if ((np.abs(zr) - rmax) <= tol).all():
            return zenv

    raise ValueError(
        "Iterative method did NOT converge."
        " You might want to increase the number of iterations and/or the tolerance."
        " The maximum slope parameter rmax is " + str(np.nanmax(np.abs(zr)))
    )

#=======================================================================================
def msg_info(message,main=False):

    if main:
       print('')
       print('='*76)
       print(' '*11 + message)
       print('='*76)
    else:
       print(message)
    print('')

