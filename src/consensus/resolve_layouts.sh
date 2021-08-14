#!/bin/bash
set -e

if [ "$#" -lt 3 ]; then
    echo "script.sh <contig layout (.paths)> <resolved mapping file> <resolved backbone file>"
    exit 239
fi

layout=$1
mapping=$2
out_path=$3
mkdir -p $(dirname $out_path)

$(dirname $(readlink -e $0))/../resolve_layouts.py --resolved-marker _i $layout $mapping > $out_path
