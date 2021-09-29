declare -i N
N=4
for LETTER in  C D E F  
do
  echo $LETTER, $N
  #sbatch spicemerge.sh $LETTER $N
  sh spicemerge.sh $LETTER $N
  N=$N-1
done

sbatch  allspicemerge.sh
