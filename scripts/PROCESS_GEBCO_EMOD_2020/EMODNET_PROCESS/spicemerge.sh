#! /bin/ksh
#SBATCH --mem=200000
#SBATCH --ntasks=1
#SBATCH --output=2output
#SBATCH --error=2error
#SBATCH --time=5

sh runmerge.sh -o /scratch/fred/REQUIRED_INPUTS -i /scratch/fred/PROCESS_EMODNET_DEC_2020/NC/ -p python

