------------------------------------------------------------------------------
Files for compiling and running AMM15 CO9_P2.0 on ARCHER2; tested October 2023
------------------------------------------------------------------------------

Based on https://code.metoffice.gov.uk/svn/roses-u/c/u/6/7/4/trunk/app/nemo_cice/rose-app.conf
! P2.0 CO9 AMM15 NEMO 404 OLD Coordinate + ERA5 + GLS; correctly ordered baltic boundary coordinates

Uses NEMO version 4.0.4 revision 13653
Runs with XIOS revision 1964

Compile using
./setup_AMM15_CO9_P2_archer -w full_path_to_work_dir -s full_path_to_P2_repository_dir

setup_P2_files.sh gives links to forcing data on ARCHER2
runamm15_872.slurm is a sample run script




