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

minimap2 -x asm5 --sam-hit-only -a -t $thread_cnt -o $name.nodes.sam $fasta $name.nodes.fasta

~/git/ngs_scripts/alignment_filter.py --min-len 50000 --min-idy 0.98 $name.nodes.sam 2> $name.filter.log | awk '{print $1}' > $name.nodes.tmp.txt
~/git/ngs_scripts/alignment_filter.py --min-len 15000 --min-idy 0.98 --query-frac .90 $name.nodes.sam 2>> $name.filter.log | awk '{print $1}' >> $name.nodes.tmp.txt
sort $name.nodes.tmp.txt | uniq > $name.nodes.txt
rm $name.nodes.tmp.txt

~/git/ngs_scripts/gfakluge/neighborhood $g $name.subgraph.gfa $name.nodes.txt 5 &> $name.nodes_rad_5.log
#done
