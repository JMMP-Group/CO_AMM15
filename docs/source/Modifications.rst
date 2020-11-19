# Modifications
================
*Last Updated: 18/11/2020*

This file describes modifications made to this NEMO configuration by way of MY_SRC
files and input datasets.

## Tidal Potential Forcing
==========================

Some aspects of the tidal potential forcing have been modified. This is based on work
done by Nicolas Bruneau. It will compile and run on NEMO versions 4.0, 4.0.1 and 4.0.2.
It has not yet been tested beyond this.

*Files Changed:*
- MY_SRC
-- tide_mod.F90
-- sbctide.F90
-- tide.h90
-- tideini.F90

*Namelist Changes:*
A new parameter `dn_love_number` is specified in the `namtide` namelist.

*Code Changes:*
> Variable love number. Previously hard coded. Set in namelist using `dn_love_number`
> Long period tidal forcing now applied. Previously set to zero despite input.
> Additional Nodal correction equation.

**Note**: Simon Mueller at NOCS worked on getting this onto the main NEMO trunk. This work
is completely but it will not be available until a future version of NEMO (perhaps 4.2).
The same applies to the improved harmonic analysis in NEMO.

## Baltic Boundary Conditions
=============================

*See: https://github.com/NOC-MSM/DEV_nibrun/wiki/AMM-series*

*Files Changed:*
- MY_SRC
-- bdydta.F90
-- bdyini.F90
-- bdy_oce.F90

*Namelist Changes:*
A new parameter `rn_ssh_shift` specified in new namelist `nambdy_ssh`.

*Description:*
Along with bathymetry changes below, the model Baltic boundary has been modified in
a number of ways to achieve better transport in the Kattegat Basin/Norwegian Trench and to
improve the stability of AMM15 to allow a larger timestep. Boundary conditions have been
changed as follows:

> BCs have been changed for the baltic boundary so that NEMO interpolates them on-the-fly.
> SSH input at boundaries can be shifted up or down using new parameter `rn_ssh_shift`.
  This happens in the three MY_SRC files listed above. A constant value is added
  to the boundary ssh in `bdydta.F90`.

## Bathymetry
==============

*See: https://github.com/NOC-MSM/NEMO_cfgs/wiki/AMM15*

*Files Changed:*

- bathy_meter.nc
- domain_cfg.nc

*Namelist Changes:* n/a

*Description:*

The bathymetry in the Kattegat region has been modified (smoothing) to allow for improved
transport in the Kattegat basin and the Norwegian Trench. Key points:

> Bathymetry has been smoothed using ROMS pythong toolbox Sikiric et al. (2009, OM). 
> Smoothing used Rmax of 0.3 for depth lower than 125m. (Kattegat)
> Baltic Sea region also smootged using Rmax = 0.1.
> Narrow channels East and West of Denmark have been enlarged, for stability.
> Manual modification to other areas, including small island near Norway and deep canyon
  near Wales.

## River Tracers
================

*See: https://github.com/NOC-MSM/NEMO_cfgs/wiki/AMM15*

*Files Changed:*

*Namelist Changes:* 
- MY_SRC
-- par_age.F90
-- trc.F90
-- trcbc.F90
-- trcini.F90
-- trcini_age.F90
-- trcnam.F90
-- trcnam_age.F90
-- trcnxt.F90
-- trcrad.F90
-- trcsbc.F90
-- trcsms_age.F90
-- trcsms_my_trc.F90
-- trcwri_age.F90

*Description:*

Modifications were made to TOP to allow for river transports and age. All modifications
are to the TOP component of NEMO. No changes are made to OCE here.

## Light Penetration
====================

*Files Changed:*
- MY_SRC
--traqsr.F90

*Namelist Changes:*

*Description:*
