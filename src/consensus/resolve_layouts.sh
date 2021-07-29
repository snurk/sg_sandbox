#!/bin/bash
set -e

if [ "$#" -lt 3 ]; then
    echo "script.sh <contig layout (.paths)> <resolved mapping file> <resolved backbone file>"
    exit 239
fi

layout=$1
mapping=$2
path=$3
mkdir -p $(dirname $path)

awk '{print $1}' $layout > $path.names
cat $layout $mapping > $path.full_mapping

$(dirname $(readlink -e $0))/../resolve_layouts.py $path.names $path.full_mapping > $path

rm $path.names
rm $path.full_mapping
