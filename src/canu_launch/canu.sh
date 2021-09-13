#!/bin/bash
set -eou

if [ "$#" -lt 3 ]; then
    echo "script.sh <out_dir> <max read len (hpc)> <reads> [additional parameters passed directly]"
    exit 239
fi

out_dir=$1
max_len=$2
reads=$3

#Options to consider
#OvlMerSize=31 (default is 22)

echo "Assembling reads $reads with Canu from $CANU_BIN into $out_dir"
echo "Additional parameters ${@:4}"

mkdir -p $out_dir
echo "$max_len" > $out_dir/max_read_len
$CANU_BIN/canu -p asm -d $out_dir -assemble canuIterationMax=1 maxInputCoverage=10000 batMemory=200 genomeSize='3.1g' gridOptionsUTGOVL='--time=12:00:00' gridOptionsmeryl='--mem-per-cpu=4g' gridOptions='--time=24:00:00' ovlThreads=32 utgOvlRefBlockLength=10000000000 ovlMerThreshold=20000 minReadLength=4000 minOverlapLength=2000 batOptions="-covgaptype uncovered -D symmetricoverlaps -stop bestedges -ef 1e-5 -readlen 4000-$max_len" stopAfter=unitig -pacbio-hifi $reads "${@:4}"
