#!/bin/bash
set -eou

if [ "$#" -lt 3 ]; then
    echo "Usage: $(basename $0) <graph.gfa> <init.gfa> <out_prefix>"
    exit 1
fi

scripts_root=$(dirname $(readlink -e $0))

graph=$1
init=$2
prefix=$3

grep "^S" $graph | cut -f 2 > $prefix.tmp
echo -e "H\tVN:Z:1.0" > $prefix.gfa
$scripts_root/../tangle-resolution/scripts/extract_subgraph.py $prefix.tmp < $init | grep "^S" | sed $'s/\tRC.*//g' >> $prefix.gfa
grep "^L" $graph >> $prefix.gfa
rm $prefix.tmp
