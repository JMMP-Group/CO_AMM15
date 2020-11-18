# Modifications
*Last Updated: 18/11/2020*

This file describes modifications made to this NEMO configuration by way of MY_SRC
files and input datasets.

## Tidal Potential Forcing

Some aspects of the tidal potential forcing have been modified. This is based on work
done by Nicolas Bruneau. It will compile and run on NEMO versions 4.0, 4.0.1 and 4.0.2.
It has not yet been tested beyond this.

*Branch:* nico_tides

*Files:*
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


## Baltic Boundary Conditions

## Bathymetry

## River Tracers