#!/bin/bash
set -e

if [ "$#" -lt 2 ]; then
    echo "script.sh <gfa> <homopolymer-compressed sequence fasta> [threads=12]"
    exit 239
fi

g=$1
fasta=$2

thread_cnt=12
if [ "$#" -ge 3 ]; then
    thread_cnt=$3
fi

#for g in $dir/simplified.gfa ; do # $dir/no_*.gfa $dir/processed.gfa ; do
name=$(basename $g .gfa)
path=$(dirname $g)/$name

if [ -f $path.nodes.fasta ]; then
    echo "File $path.nodes.fasta was found in the folder and will be reused"
    ln -s $path.nodes.fasta
else
    awk '/^S/{print ">"$2"\n"$3}' $g | fold > $name.nodes.fasta
fi

minimap2 -I 8G -x asm5 -a -t $thread_cnt -o $name.nodes.sam $name.nodes.fasta $fasta

#~/git/ngs_scripts/alignment_filter.py --min-len 50000 --min-idy 0.98 $name.nodes.sam 2> $name.filter.log | awk '{print $1}' > $name.nodes.tmp.txt
~/git/ngs_scripts/alignment_filter.py --min-idy 0.99 --query-frac .50 $name.nodes.sam 2>> $name.filter.log > $name.align.info

#filtering out queries which lie within the unitig
awk '{if ($3 > 0.99) print $1}' $name.align.info | sort > $name.contained_queries.txt

#NB only works with query names that can not clash with anything else
grep -v -f $name.contained_queries.txt $name.align.info | awk '{print $5}' | sort | uniq > $name.nodes.txt

~/git/ngs_scripts/gfakluge/neighborhood $g $name.subgraph.gfa $name.nodes.txt 5 &> $name.nodes_rad_5.log
