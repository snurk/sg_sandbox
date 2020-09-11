#!/bin/bash
set -e

if [ "$#" -lt 2 ]; then
    echo "script.sh <subgraph gfa> <homopolymer-compressed reads> [threads=12]"
    exit 239
fi

subgraph=$1
name=$(dirname $1)/$(basename $1 .gfa)
reads=$2

thread_cnt=12
if [ "$#" -ge 3 ]; then
    thread_cnt=$3
fi

if [ -f ${name}.nodes.fasta ];
then
    echo "Node sequences already extracted"
else
    echo "Extracting node sequences"
    awk '/^S/{print ">"$2"\n"$3}' $subgraph | fold > ${name}.nodes.fasta
fi

echo "Aligning reads with minimap2"
#This will exclude unmapped (0x4), secondary (0x100), and supplementary (0x800) alignments.
minimap2 -x asm5 -a -I 10G -t $((thread_cnt - 1)) ${name}.nodes.fasta $reads 2> ${name}.read_align.log | samtools view -F 2308 -q 40 -b > ${name}.reads.bam

echo "Filtering read alignments"
~/git/ngs_scripts/alignment_filter.py --min-len 5000 --query-frac 0.9 --min-idy 0.99 --filtered ${name}.reads.filtered.bam ${name}.reads.bam &> ${name}.read_filter.log

#echo "Sorting read alignments"
#samtools sort ${name}.reads.filtered.bam --threads $((thread_cnt - 1)) -o ${name}.reads.sorted.bam

echo "Computing coverage"
~/git/ngs_scripts/alignment_stats/compute_coverage.py ${name}.reads.filtered.bam &> ${name}.cov
