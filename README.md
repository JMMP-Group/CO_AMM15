# CO9_AMM15
CO9 configuration of the Atlantic Margin Model (1.5 km resolution)

This repository is still under construction. It has begun to be populated so that version tracking and branching can begin. 

Directories:

	- EXP00
		Main, combined run directory
	- EXP_nico
		Run directory used by Nico for 10-year RIVER run

	- MY_SRC
		Main, combined MY_SRC directory
	- MY_SRC_river15
		MY_SRC for freshwater tracers from river runoff, Baltic (Nico)
	- MY_SRC_tides
		MY_SRC for modifications to tide forcing.
		
	- arch
		ARCHER architecture files for NEMO and XIOS
	- scripts
		Miscellaneous scripts for running, building and moving configurations/files

	- cpp_co9-amm15.fcm
		NEMO compile/build keys.

Input Files:

These can be found on a few different locations in /projectsa/, ARCHER and JASMIN:

	- Domain
	- Surface Forcing
	- Lateral Boundary Forcing
	- River Forcing
	- Other
