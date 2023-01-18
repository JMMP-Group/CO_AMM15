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

Firther background to CO9-AMM15 can be found on the wiki:

https://github.com/JMMP-Group/CO_AMM15/wiki

## Configuration Input Files

|  **Input** | **Download Location** |
|-------------- | -------------- |
| **Final Bathy used to make below domain file** | https://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/AMM15_BATHY/G-E-G-NICO_BLOCK_River_10mMIN.nc |
| **P1.5**  **ME** **Domain_cfg.nc** | http://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/DOMAIN_CFG/domain_cfg_sig9_itr3_MEs_01-002_350-1400_local_opt_v1.nc	 |
| **P1.5b** **SF12** **Domain_cfg.nc** | http://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/DOMAIN_CFG/GEG_SF12.nc	 |
| **AMM15 Coordinates files** | [http://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/COORDINATES/amm15.coordinates.rim15.nc](https://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/COORDINATES/amm15.bdy.coordinates.rim15.nc)	 |
| **Open ocean boundary coordinates.bdy.nc** | [http://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/COORDINATES/amm15.bdy.coordinates.rim15.nc](https://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/COORDINATES/amm15.bdy.coordinates.rim15.nc)	 |
| **Baltic coordinates.bdy.nc** | https://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/COORDINATES/amm15.baltic.bdy.coordinates.nc |
---

## Sample Forcing Files

| **Forcing** | **Download Location** |
|-------------- | ------------------|
| **Surface boundary** |[ http://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/FORCING/SBC/ERA5/](https://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/FORCING/SBC/ERA5/) |
| **Open ocean boundary** |[ http://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/FORCING/BDY/](https://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/FORCING/BDY/) |
| **Atlantic Baroclinic No vertical intrpolation** | [http://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/FORCING/BDY/EXPER_NO_VERT_BDY_SJPZ_A_AND_D/](https://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/FORCING/BDY/EXPER_NO_VERT_BDY_SJPZ_A_AND_D/) |
| **Atlantic Barotropic** | https://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/FORCING/BDY/SJPZ_A_AND_D_BT/ |
| **Baltic boundary** | [https://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/FORCING/BDY/amm15_Baltic/ |
| **River runoff** | http://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/FORCING/RIVERS/ |
| **Tide** | [http://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/FORCING/TIDES/FES2014/](https://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/FORCING/TIDES/FES2014/) |
| **Initial condition** | http://gws-access.jasmin.ac.uk/public/jmmp_collab/AMM15/inputs/IC/ |
