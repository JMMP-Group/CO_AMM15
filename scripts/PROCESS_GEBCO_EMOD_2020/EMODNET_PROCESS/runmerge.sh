
#bash script

usage() { echo "Usage: $0 [-i <input_directory>] [-o <output_directory>] [-p <python command>]" 1>&2; exit 1; }

while getopts ":i:o:p:" t; do
    case "${t}" in
        i)
            i=${OPTARG}
            ;;
        o)
            o=${OPTARG}
            ;;
        p)
            p=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${o}" ] ; then
     echo "\n-o <output_dir> not specified\n"
     usage
fi
if [ -z "${i}" ]; then
     echo "\n-i <input_dir> not specified\n"
     usage
fi
if [ -z "${p}" ]; then
     echo "\n-o <python command> not specified\n"
     usage
fi

if [ ! -d "${o}" ] ; then
     echo "\n ${o} <output_dir> does not exist\n"
     usage
fi
if [ ! -d "${i}" ]; then
     echo "\n ${i} <input_dir> does not exist\n"
     usage
fi

echo "<input_dir> = ${i}"
echo "<output_dir> = ${o}"
echo "python = ${p}"



declare -i N
N=4
for LETTER in  C D E F  
do
  echo $LETTER, $N
  # merge along lon
  ${p} merge_xarray.py    ${i}/$LETTER ${o}/$N
  N=($N-1)
done

# merge along lat
${p}   final_merge.py ${o}
