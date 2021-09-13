#!/bin/bash
set -eou

if [ "$#" -lt 2 ]; then
    echo "script.sh <assembly folder> <output .fasta[.gz]>"
    exit 239
fi

assembly=$1
out=$2

$CANU_BIN/fixErrors -S $assembly/asm.seqStore/ -red $assembly/unitigging/3-overlapErrorAdjustment/red.red -O $out
