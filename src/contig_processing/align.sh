#!/bin/bash
set -e

if [ "$#" -lt 2 ]; then
    echo "script.sh <query> <ref> <bam> [thread count = 8]"
    exit
fi

query=$1
ref=$2
bam=$3

threads=8
if [ "$#" -ge 4 ]; then
    threads=$4
fi

module load minimap2
module load samtools

minimap2 -a -H -t $threads -x asm20 -I 8G $ref $query | samtools view -b -F 4 > $bam
