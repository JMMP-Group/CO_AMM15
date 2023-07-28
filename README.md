# CO_AMM15

Coastal Ocean (CO) configuration of the Atlantic Margin Model (1.5 km resolution)

[Citations for releases can be found on Zenodo](https://doi.org/10.5281/zenodo.6497328)

The Joint Marine Modelling Programme [(JMMP)](https://www.metoffice.gov.uk/research/approach/collaboration/joint-marine-modelling-programme) provides world-class and easily accessible national capability, ocean modelling infrastructure and configurations to support the UK’s scientific research and operational prediction systems for ocean, weather and climate. It is partnership between the Met Office and British Antarctic Survey, National Oceanography Centre and Centre for Polar Observation and Modelling.

Model configurations are underpinned by the Nucleus for European Modelling of the Ocean [(NEMO)](https://www.nemo-ocean.eu) framework. JMMP works closely with the NEMO consortium to develop the underpinning model capability. 

---

At present, the main is for CO10 development. CO9 will continute to progress under the CO9 branch, with associated v9.x.x releases. See tags for the latest release.

Directories:

    - ARCHER2
        ARCHER2 specific EXP files
    - BUILD_CFG
        scrips and files related to specific checkpoints
	- EXP00
		Main experiment
	- MY_SRC
		Custom NEMO 4.0.4 scripts
	- arch
		ARCHER architecture files for NEMO and XIOS
    - docs
        Miscellaneous documentation
	- scripts
		Miscellaneous scripts for running, building and moving configurations/files
	- cpp_co9-amm15.fcm
		NEMO compile/build keys.

Further background to CO-AMM15 can be found on the wiki:

https://github.com/JMMP-Group/CO_AMM15/wiki

## Configuration Input Files

|  **Input** | **Download Location** |
|-------------- | -------------- |
| **Final Bathy used to make below domain file** | https://gws-access.jasmin.ac.uk/public/jmmp/AMM15/AMM15_BATHY/G-E-G-NICO_BLOCK_River_10mMIN.nc |
| **P1.5**  **ME** **Domain_cfg.nc** | http://gws-access.jasmin.ac.uk/public/jmmp/AMM15/DOMAIN_CFG/domain_cfg_sig9_itr3_MEs_01-002_350-1400_local_opt_v1.nc	 |
| **P1.5b** **SF12** **Domain_cfg.nc** | http://gws-access.jasmin.ac.uk/public/jmmp/AMM15/DOMAIN_CFG/GEG_SF12.nc	 |
| **AMM15 Coordinates files** | [http://gws-access.jasmin.ac.uk/public/jmmp/AMM15/COORDINATES/amm15.coordinates.nc](https://gws-access.jasmin.ac.uk/public/jmmp/AMM15/COORDINATES/amm15.coordinates.nc)	 |
| **Open ocean boundary coordinates.bdy.nc** | [http://gws-access.jasmin.ac.uk/public/jmmp/AMM15/COORDINATES/amm15.bdy.coordinates.rim15.nc](https://gws-access.jasmin.ac.uk/public/jmmp/AMM15/COORDINATES/amm15.bdy.coordinates.rim15.nc)	 |
| **Baltic coordinates.bdy.nc** | https://gws-access.jasmin.ac.uk/public/jmmp/AMM15/COORDINATES/amm15.baltic.bdy.coordinates.nc |
| **XIOS xml files** | https://gws-access.jasmin.ac.uk/public/jmmp/AMM15/XML/ |
---

## Sample Forcing Files

| **Forcing** | **Download Location** |
|-------------- | ------------------|
| **Surface boundary** |[ http://gws-access.jasmin.ac.uk/public/jmmp/AMM15/FORCING/SBC/ERA5/](https://gws-access.jasmin.ac.uk/public/jmmp/AMM15/FORCING/SBC/ERA5/) |
| **Open ocean boundary** |[ http://gws-access.jasmin.ac.uk/public/jmmp/AMM15/FORCING/BDY/](https://gws-access.jasmin.ac.uk/public/jmmp/AMM15/FORCING/BDY/) |
| **Atlantic Baroclinic No vertical intrpolation** | [http://gws-access.jasmin.ac.uk/public/jmmp/AMM15/FORCING/BDY/EXPER_NO_VERT_BDY_SJPZ_A_AND_D/](https://gws-access.jasmin.ac.uk/public/jmmp/AMM15/FORCING/BDY/EXPER_NO_VERT_BDY_SJPZ_A_AND_D/) |
| **Atlantic Barotropic** | https://gws-access.jasmin.ac.uk/public/jmmp/AMM15/FORCING/BDY/SJPZ_A_AND_D_BT/ |
| **Baltic boundary** | [https://gws-access.jasmin.ac.uk/public/jmmp/AMM15/FORCING/BDY/amm15_Baltic/ |
| **River runoff** | http://gws-access.jasmin.ac.uk/public/jmmp/AMM15/FORCING/RIVERS/ |
| **Tide** | [http://gws-access.jasmin.ac.uk/public/jmmp/AMM15/FORCING/TIDES/FES2014/](https://gws-access.jasmin.ac.uk/public/jmmp/AMM15/FORCING/TIDES/FES2014/) |
| **ME Initial condition** (P1.5) | http://gws-access.jasmin.ac.uk/public/jmmp/AMM15/https://gws-access.jasmin.ac.uk/public/jmmp/AMM15/RESTARTS/RESTART_BASED_ONCO7_TO_GEG_NICO_BALTIC_BLOCK_BUT_10M_MIN_RIV_DEP_MEs_01-002_350-1400/
| **SF12 Initial Condition (P1.5b)** | https://gws-access.jasmin.ac.uk/public/jmmp/AMM15/RESTARTS/RESTART_BASED_ONCO7_20040101_TO_GEG_NICO_BALTIC_BLOCK_BUT_10M_MIN_RIV_DEP/ |
