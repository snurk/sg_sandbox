#!/bin/bash
set -eou

if [ "$#" -lt 3 ]; then
    echo "script.sh <canu_exec> <out_dir> <reads> [additional parameters passed directly]"
    exit 239
fi

canu_exec=$1
out_dir=$2
max_len=$3
reads=$4

#Options to consider
#OvlMerSize=31 (default is 22)

echo "Assembling reads $reads with Canu from $canu_exec into $out_dir"
echo "Additional parameters ${@:4}"

$canu_exec -p asm -d $out_dir -assemble canuIterationMax=1 maxInputCoverage=10000 batMemory=200 genomeSize='3.1g' gridOptionsUTGOVL='--time=12:00:00' gridOptionsmeryl='--mem-per-cpu=4g' gridOptions='--time=24:00:00' ovlThreads=32 utgOvlRefBlockLength=10000000000 ovlMerThreshold=20000 minReadLength=4000 minOverlapLength=2000 batOptions="-covgaptype uncovered -D symmetricoverlaps -stop bestedges -ef 1e-5 -readlen 4000-$max_len" stopAfter=unitig -pacbio-hifi $reads "${@:4}"
