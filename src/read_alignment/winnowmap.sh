#!/bin/bash
set -e

if [ "$#" -lt 3 ]; then
    echo "script.sh <query> <ref> <bad_kmers> <bam> [thread count = 8]"
    exit
fi

query=$1
ref=$2
bad_kmers=$3
bam=$4

threads=8
if [ "$#" -ge 5 ]; then
    threads=$5
fi

module load samtools

winnowmap -W $bad_kmers -a -x asm -I 8G -t $threads $ref $query | samtools view -b -F 4 > $bam
