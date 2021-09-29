#!/bin/bash



usage() { echo "Usage: $0 [-i <input_directory>] [-o <output_directory>]" 1>&2; exit 1; }

while getopts ":i:o:" p; do
    case "${p}" in
        i)
            i=${OPTARG}
            ;;
        o)
            o=${OPTARG}
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
ZIPFILES=`ls ${i}/*.zip`
for file in  $ZIPFILES
do
	ls $file
	unzip -o $file -d ${o}  &
done



