#!/bin/bash
set -eou

if [ "$#" -lt 2 ]; then
    echo "script.sh <out_dir> <reads> [additional parameters passed directly]"
    exit 239
fi

out_dir=$1
reads=$2

#Options to consider
#OvlMerSize=31 (default is 22)

echo "Assembling reads $reads with Canu from $CANU_BIN into $out_dir"
echo "Bogart will mark reads longer than $MAX_BOGART_LEN in hpc as 'bad'"
echo "Additional parameters ${@:3}"

mkdir -p $out_dir
$CANU_BIN/canu -p asm -d $out_dir -assemble canuIterationMax=1 maxInputCoverage=10000 batMemory=200 genomeSize='3.1g' gridOptionsUTGOVL='--time=12:00:00' gridOptionsmeryl='--mem-per-cpu=4g' gridOptions='--time=24:00:00' ovlThreads=32 utgOvlRefBlockLength=10000000000 ovlMerThreshold=20000 minReadLength=4000 minOverlapLength=1000 batOptions="-covgaptype uncovered -D symmetricoverlaps -stop bestedges -ef 1e-5 -readlen 1-$MAX_BOGART_LEN" stopAfter=unitig -pacbio-hifi $reads "${@:3}"

#For debug using E.coli
#$CANU_BIN/canu -p asm -d $out_dir -assemble canuIterationMax=1 maxInputCoverage=10000 batMemory=20 genomeSize='4.5m' gridOptionsUTGOVL='--time=12:00:00' gridOptionsmeryl='--mem-per-cpu=4g' gridOptions='--time=4:00:00' ovlThreads=16 utgOvlRefBlockLength=10000000000 ovlMerThreshold=20000 minReadLength=4000 minOverlapLength=1000 batOptions="-covgaptype uncovered -D symmetricoverlaps -stop bestedges -ef 1e-5 -readlen 1-$MAX_BOGART_LEN" stopAfter=unitig -pacbio-hifi $reads "${@:3}"
#$CANU_BIN/canu -p asm -d $out_dir -assemble useGrid=false canuIterationMax=1 maxInputCoverage=10000 batMemory=20 genomeSize='4.5m' ovlThreads=1 corThreads=1 redThreads=1 oeaThreads=1 batThreads=1 minReadLength=4000 minOverlapLength=1000 batOptions="-covgaptype uncovered -D symmetricoverlaps -stop bestedges -ef 1e-5 -readlen 1-$MAX_BOGART_LEN" stopAfter=unitig -pacbio-hifi $reads "${@:3}"
