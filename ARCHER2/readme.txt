Files for compiling and running AMM15 on ARCHER2.

Module list (current at 18 July 2022)

module swap PrgEnv-cray/8.0.0 PrgEnv-gnu/8.1.0
module swap craype-network-ofi craype-network-ucx
module swap cray-mpich cray-mpich-ucx
module load cray-hdf5-parallel/1.12.0.7
module load cray-netcdf-hdf5parallel/4.7.4.7
module load libfabric

Based on https://code.metoffice.gov.uk/svn/roses-u/c/d/6/1/1/trunk/app/nemo_cice/rose-app.conf
! P0.9b CO9 AMM15 NEMO 404 ERA5+update GLS options + Atlantic 3D UV specified,Tra specified,LN_FULL_TRUE +TPXO (u-cd513) ; 

Uses NEMO version NERC/NEMO_4.0.4_CO9_package_tides@r14364
Runs with XIOS revision 1903

May 2023
Added setup script for compiling AMM15_CO9_1.5 on ARCHER2: setup_AMM15_CO9_1.5

