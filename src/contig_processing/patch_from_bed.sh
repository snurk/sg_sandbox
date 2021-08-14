#!/bin/bash
set -e

if [ "$#" -lt 3 ]; then
    echo "script.sh <fragments.fasta> <patch.bed> <out_name>"
    exit 239
fi

name=$3

echo ">$name" > tmp.fasta
bedtools getfasta -s -fi $1 -bed $2 | grep -v ">" >> tmp.fasta

$(dirname $(readlink -e $0))/contig_length_filter.py 1 tmp.fasta $name.fasta

rm -f tmp.fasta
