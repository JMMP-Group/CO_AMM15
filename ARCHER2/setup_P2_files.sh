# link to input files (other filenames are specified in namelist_cfg)
# bdy coordinate files
ln -s /work/n01/n01/shared/CO_AMM15/P1_INPUTS/COORDINATES/amm15.bdy.coordinates.rim15.nc coordinates.bdy.nc 
ln -s /work/n01/n01/shared/CO_AMM15/P1_INPUTS/COORDINATES/amm15_baltic.bdy9.coordinates.nc coordinates.skagbdy.nc

# RESTART FOR P1.5b
ln -s /work/n01/n01/shared/CO_AMM15/P1_INPUTS/FORCING/RESTART/RESTART_BASED_ONCO7_20040101_TO_GEG_NICO_BALTIC_BLOCK_BUT_10M_MIN_RIV_DEP INITIAL_RESTART

#RIVERS
ln -s /work/n01/n01/shared/CO_AMM15/P1_INPUTS/FORCING/RIVERS/Spread_Loirre_Garron_N_G.nc rivers.nc

# BDY, TIDE and ERA5 forcing
ln -s /work/n01/n01/shared/CO_AMM15/P1_INPUTS/FORCING .

# make sure directory for outputting RESTART files exists
mkdir RESTART
