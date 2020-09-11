#!/bin/bash
set -e

if [ "$#" -lt 2 ]; then
    echo "script.sh <contig layout> <resolved mapping file>"
    exit 239
fi

layout=$1

awk '{print $1}' $layout > contig_names.txt
cat $layout $2 > full_mapping.txt

$(dirname $0)/resolve_layouts.py contig_names.txt full_mapping.txt > backbone_layout.txt
