#! /bin/ksh
#SBATCH --mem=100000
#SBATCH --ntasks=1
#SBATCH --output=expand_output
#SBATCH --error=expand_error
#SBATCH --time=15

echo "python EXPAND_AMM15_CUBE.py  -a $1 -i $2  -o  $3"
python EXPAND_AMM15_CUBE.py  -a $1 -i $2  -o  $3

