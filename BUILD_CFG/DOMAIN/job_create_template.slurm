#!/bin/bash
#SBATCH --job-name=make_domaincfg
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=120
#SBATCH --cpus-per-task=1
#SBATCH --account=n01-ACCORD
#SBATCH --partition=standard
#SBATCH --reservation=shortqos
#SBATCH --qos=short

# Setup the job environment (this module needs to be loaded before any other modules)
module -s restore /work/n01/shared/acc/n01_modules/ucx_env

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1
#export UCX_IB_REG_METHODS=direct

# Launch the parallel job
#   srun picks up the distribution from the sbatch options
srun --distribution=block:block --hint=nomultithread XXX_TDIR_XXX/DOMAINcfg/BLD/bin/make_domain_cfg.exe

