#!/bin/bash
set -eou

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <reference> <chr or contig name to select> <list of reads>"
    exit 1
fi

base_path=$(dirname $(readlink -e $0))

ref=$1
chr=$2

>input.fofn
count=0
for i in `seq 3 $#`; do
   echo "${!i}" >> input.fofn
   ((count=count+1))
done
out=$(basename $ref .fasta)"_$chr"

mkdir -p $out/logs
mkdir -p $out/map

if [ ! -f $out/merged.sorted.bam ] ; then
   sbatch -o $out/logs/map.%A_%a.out -D `pwd` --mem=100G --cpus-per-task=40 --time=36:00:00 -a 1-$count $base_path/winnowmap_ont.sh $ref $out/map > map.jid
   
   WAIT="afterok:"`cat map.jid | tr '\n' ',afterok:'`
   WAIT=${WAIT%,}
   sbatch -o $out/logs/subset.%A_%a.out -D `pwd` --mem=100G --cpus-per-task=8 --time=36:00:00 --dependency=$WAIT $base_path/subset_winnowmap_ont.sh $chr $out 
fi
