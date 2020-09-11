#!/bin/bash
set -e

if [ "$#" -lt 4 ]; then
    echo "script.sh <chr> <contig info file> <reference alignment tab file> <out bed>"
    exit 239
fi

chr=$1
ctg_info=$2
ref_align=$3
out_bed=$4

#tmp file is appended to deal with insertion of gap in acrocentrics
$(dirname $0)/gap_patch.py $ctg_info $ref_align | sed 's/DJ/a/g' | sed 's/PJ/b/g' >> $out_bed.tmp

sort $out_bed.tmp | uniq | sed 's/.*: //g' | sed 's/\.a/\.DJ/g' | sed 's/\.b/\.PJ/g' > $out_bed
#rm -f $out_bed.tmp
