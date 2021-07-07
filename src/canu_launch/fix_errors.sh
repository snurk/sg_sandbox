#!/bin/bash
set -eou

if [ "$#" -lt 3 ]; then
    echo "script.sh <fix_exec> <assembly folder> <output .fasta[.gz]>"
    exit 239
fi

fix_exec=$1
assembly=$2
out=$3

$1 -S $assembly/asm.seqStore/ -red $assembly/unitigging/3-overlapErrorAdjustment/red.red -O $out
