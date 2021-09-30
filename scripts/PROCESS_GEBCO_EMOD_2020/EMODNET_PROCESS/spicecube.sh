#! /bin/ksh
#SBATCH --mem=200000
#SBATCH --ntasks=1
#SBATCH --output=2output
#SBATCH --error=2error
#SBATCH --time=5

python MAKE_EMODNET_CUBE.py -i $1 -o  $2

