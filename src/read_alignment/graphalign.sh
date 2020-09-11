#!/bin/bash
set -e

if [ "$#" -lt 3 ]; then
    echo "script.sh <gfa> <in fasta> <out (.gaf/.gam/.json)> [threads=16]"
    exit 239
fi

gfa=$1
in_fasta=$2
out=$3

#chrY.simplified.compressed.gfa
thread_cnt=16
if [ "$#" -ge 4 ]; then
    thread_cnt=$4
fi

~/git/GraphAligner/bin/GraphAligner -t $thread_cnt -x dbg -g $gfa -f $in_fasta -a tmp.$out
mv tmp.$out $out
