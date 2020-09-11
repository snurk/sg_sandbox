#!/bin/bash
set -e

if [ "$#" -lt 1 ]; then
    echo "script.sh <chromosome name>"
    echo "telomere_patch.bed files should be in the folder"
    exit 239
fi

chr=$1

echo ">$chr" > noformat.fasta
bedtools getfasta -s -fi ../for_patch.fasta -bed telomere_patch.bed | grep -v ">" >> noformat.fasta

~/git/ngs_scripts/contig_processing/contig_length_filter.py 1 noformat.fasta $chr.fasta
