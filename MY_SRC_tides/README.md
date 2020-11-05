# MY_SRC

Contains developement code for SE-NEMO
Square brackets show comment markers in the code to find and identify changes.

Look for initials: NB, davbyr

sbctide.F90
   - [NB] Applies Love Number from namelist.
   - [NB] Applies long period tide forcing (previously set to zero).

tide.h90
   - [NB] Updated tidal potential forcing data. More constituents, Schureman method.

tide_mod.F90
   - [NB] Updated nodal factor equation.

tideini.F90
   - [NB] Reads Love number from namelist and outputs to Ocean.output.
