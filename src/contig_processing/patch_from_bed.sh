#!/bin/bash
set -eou

if [ "$#" -lt 3 ]; then
    echo "script.sh <fragments.fasta> <patch.bed> <prefix>"
    exit 239
fi

prefix=$3
name=$(basename $prefix)

echo ">$name" > $prefix.tmp
bedtools getfasta -s -fi $1 -bed $2 | grep -v ">" >> $prefix.tmp

$(dirname $(readlink -e $0))/length_filter.py $prefix.tmp > $prefix.fasta

rm -f $prefix.tmp
