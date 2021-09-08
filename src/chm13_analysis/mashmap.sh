#!/bin/bash
set -eou

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <query> <reference> <out> [identity threshold, default = 95] [threads, default = 16]"
    exit 1
fi

query=$1
reference=$2
out=$3

identity=95
if [ "$#" -gt 3 ]; then
    identity=$4
fi

threads=16
if [ "$#" -gt 4 ]; then
    threads=$5
fi

module load mashmap

mashmap -f map -s 10000 -t $threads --pi $identity -q $query -r $reference -o $out
