#!/bin/bash
set -eou

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <graph> <reads> <prefix> [threads, default = 8]"
    exit 1
fi

graph=$1
reads=$2
prefix=$3

threads=8
if [ "$#" -gt 3 ]; then
    threads=$4
fi

echo "Aligning reads $reads to graph $graph (in $threads threads). Output in $prefix.gaf, log in $prefix.ga_align.log"

/home/nurks2/git/GraphAligner/bin/GraphAligner -x vg -b 50 --multiseed-DP 1 --X-drop 1000000 --precise-clipping 0.9 --multimap-score-fraction 1 --min-alignment-score 10000 -t $threads -g $graph -f $reads -a $prefix.tmp.gaf &> $prefix.ga_align.log
mv $prefix.tmp.gaf $prefix.gaf
