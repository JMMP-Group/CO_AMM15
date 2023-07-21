Files for compiling and running AMM15 on ARCHER2.

Module list (current at 21 July 2023)

GNU 
module swap PrgEnv-cray PrgEnv-gnu
module swap craype-network-ofi craype-network-ucx
module swap cray-mpich cray-mpich-ucx
module load cray-hdf5-parallel/1.12.2.1
module load cray-netcdf-hdf5parallel/4.9.0.1
module load libfabric

Based on https://code.metoffice.gov.uk/svn/roses-u/c/d/6/1/1/trunk/app/nemo_cice/rose-app.conf
! P0.9b CO9 AMM15 NEMO 404 ERA5+update GLS options + Atlantic 3D UV specified,Tra specified,LN_FULL_TRUE +TPXO (u-cd513) ; 

Uses NEMO version NERC/NEMO_4.0.4_CO9_package_tides@r14364
Runs with XIOS revision 1903

May 2023
Added setup script for compiling AMM15_CO9_1.5 on ARCHER2: setup_AMM15_CO9_1.5

