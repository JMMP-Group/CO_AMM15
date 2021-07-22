#!/bin/bash

################################
################################
# make_domain_cfg.sh
################################
################################
#
# This script is generates s-sigma vertical coordinates with
#  the provided coordinates and bathymetry netCDF files.

## Make paths
################################
export WORK=/work/n01/n01/$USER
export WDIR=$WORK/CO9_AMM15 # Also the git clone directory
export NEMO_VER=4.0.4
export NEMO=$WDIR/BUILD_CFG/$NEMO_VER
#export CDIR=$NEMO/cfgs
#export EXP=$WDIR/RUN_DIRECTORIES/EXPconstTS
#export XIOS_DIR=$WDIR/BUILD_EXE/XIOS/xios-2.5
#export XIOS1_DIR=$WDIR/BUILD_EXE/XIOS/xios-1.0

export TDIR=$NEMO/tools
export DOMAIN=$WDIR/BUILD_CFG/DOMAIN
export DOWNLOADS=$WDIR/DOWNLOADS



if [ ! -d "$DOMAIN" ]; then
  mkdir $WDIR/BUILD_CFG
  mkdir $WDIR/BUILD_CFG/DOMAIN
fi

if [ ! -d "$NEMO_VER" ]; then
  mkdir $NEMO
fi

if [ ! -d "$TDIR" ]; then
  mkdir $TDIR
fi

if [ ! -d "$DOWNLOADS" ]; then
  mkdir $WDIR/DOWNLOADS
fi


# Build tools
################################
# This is a bit combersome and I am not svn savvy os download NEMO code in order to get the tools working

echo "Checking out NEMO repository"
case "${NEMO_VER}"
  in
  4.0.4)   echo "NEMO Verion 4.0.6 will be checked out"
           ;;
  *)       echo "NEMO Version not recognised"
           echo "Versions available at present: 4.0.6"
           exit 1
esac
svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r$NEMO_VER --depth empty $NEMO
svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r$NEMO_VER/src --depth infinity $NEMO/src
svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r$NEMO_VER/cfgs/SHARED $NEMO/cfgs/SHARED
svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r$NEMO_VER/cfgs/AMM12 $NEMO/cfgs/AMM12
svn export http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r$NEMO_VER/cfgs/ref_cfgs.txt $NEMO/cfgs/ref_cfgs.txt

cd $NEMO
# Now check EXTERNALS revision number before checking out the rest
for ext_name in mk FCM IOIPSL
  do
  ext=`svn propget svn:externals | grep $ext_name | cut -c2-`
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/$ext
done

ext=`svn propget svn:externals | grep makenemo | cut -c2-`
svn export http://forge.ipsl.jussieu.fr/nemo/svn/$ext

cd $NEMO/ext/FCM/lib/Fcm
sed -e "s/FC_MODSEARCH => '',  /FC_MODSEARCH => '-J',/" Config.pm > tmp_file
mv tmp_file Config.pm

# copy the appropriate architecture file into place
mkdir $NEMO/arch
cp $WDIR/arch/nemo/arch-X86_ARCHER2-Cray.fcm $NEMO/arch/arch-X86_ARCHER2-Cray.fcm

## Edit ARCH file - not required for this tools only build
## Dirty fix to hard wire path otherwise user will have to set XIOS_DIR in every new shell session
#sed "s?XXX_XIOS_DIR_XXX?$XIOS_DIR?" $NEMO/arch/arch-X86_ARCHER2-Cray.fcm > tmp_arch
#mv tmp_arch $NEMO/arch/arch-X86_ARCHER2-Cray.fcm

## Now we actually get the tools
cd $NEMO	
svn co http://forge.ipsl.jussieu.fr/nemo/svn/NEMO/branches/UKMO/tools_r4.0-HEAD_dev_MEs $NEMO/tools


#load modules
module -s restore /work/n01/shared/acc/n01_modules/ucx_env

# compile tools
cd $TDIR
./maketools -m X86_ARCHER2-Cray -n REBUILD_NEMO
./maketools -m X86_ARCHER2-Cray -n DOMAINcfg


# Execute tool
################################

  ## Obtain the appropriate namelist (modify it if necessary)
  # Multi-Envelope vertical coordinates
  cp $DOMAIN/ME_DOMAINcfg_namelist_cfg $TDIR/DOMAINcfg/namelist_cfg

  # Ensure that the namelist_cfg has the appropriate parameters and number of lat,lon,depth levels set

  # Ensure the coordinates and bathymetry files, previously generated, are in place.
  ln -s $DOWNLOADS/amm15.coordinates.nc $TDIR/DOMAINcfg/coordinates.nc
  ln -s $DOWNLOADS/amm15.10m.bathy_meter.nc $TDIR/DOMAINcfg/bathy_meter.nc

  ## Make an adjustment to the DOMAINcfg source code to accomodate more varied vertical coords.
  ## Done in make_tools.sh
  #cp $DOMAIN/domzgr.f90.melange $TDIR/DOMAINcfg/src/domzgr.f90

  # Edit job script
  sed "s?XXX_TDIR_XXX?$TDIR?" $DOMAIN/job_create_template.slurm > $TDIR/DOMAINcfg/job_create.slurm
  
  # Submit the domain creation as a job,
  cd $TDIR/DOMAINcfg
  sbatch job_create.slurm
  
  #wait for domain creation job to finish
  for i in {0..7}; do #8 tiles
  while [ ! -f domain_cfg_000$i.nc ] ;
  do
      echo  "wait for domain creation job to finish"
      sleep 60
  done
  done
  
  # Rebuild the files. Here there are 8 tiles (and rebuilding on a single thread) 
  $TDIR/REBUILD_NEMO/rebuild_nemo -t 1 domain_cfg 8

  # After create copy it and store it for further use
  cp $TDIR/DOMAINcfg/domain_cfg.nc $DOMAIN/domain_cfg_AMM15.nc
  rm $TDIR/DOMAINcfg/domain_cfg_000*.nc #remove tiles
  cd $WORK
