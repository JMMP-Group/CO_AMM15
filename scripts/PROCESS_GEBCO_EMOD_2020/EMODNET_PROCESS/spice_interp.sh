#! /bin/ksh
#SBATCH --mem=100000
#SBATCH --ntasks=1
#SBATCH --output=2output
#SBATCH --error=2error
#SBATCH --time=15


python EXPAND_AMM15_CUBE.py -a $1  -c $2 -o $3
