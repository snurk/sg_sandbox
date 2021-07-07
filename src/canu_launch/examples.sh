#!/bin/bash
set -eou

base_path=$(dirname $(readlink -e $0))
canu_path=~/git/canu/build/bin
export CANU_BIN=$(readlink -e $canu_path)

#example 1 -- single read file
$base_path/canu.sh $CANU_BIN/canu out_canu reads.fa.gz onSuccess="$base_path/on_init_complete.sh"

#example 2 -- multiple read files
$base_path/canu.sh $CANU_BIN/canu out_canu reads1.fastq.gz -pacbio-hifi reads2.fastq.gz -pacbio-hifi reads3.fastq.gz onSuccess="$base_path/on_init_complete.sh"
