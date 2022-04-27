# CO_AMM15
Coastal Ocean (CO) configuration of the Atlantic Margin Model (1.5 km resolution)

The Joint Marine Modelling Programme [(JMMP)](https://www.metoffice.gov.uk/research/approach/collaboration/joint-marine-modelling-programme) provides world-class and easily accessible national capability, ocean modelling infrastructure and configurations to support the UK’s scientific research and operational prediction systems for ocean, weather and climate. It is partnership between the Met Office and British Antarctic Survey, National Oceanography Centre and Centre for Polar Observation and Modelling.

Model configurations are underpinned by the Nucleus for European Modelling of the Ocean [(NEMO)](https://www.nemo-ocean.eu) framework. JMMP works closely with the NEMO consortium to develop the underpinning model capability. 

---

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
