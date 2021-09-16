#!/bin/bash
set -eou

base_path=$(dirname $(readlink -e $0))
#/home/nurks2/git/canu/build/bin
export CANU_BIN=$(readlink -e $base_path/../../canu/build/bin/)
export MAX_HPC_LEN=18000

#example 1 -- single read file
$base_path/canu.sh out_canu reads.fa.gz onSuccess="$base_path/on_init_complete.sh"

#example 2 -- multiple read files
$base_path/canu.sh out_canu reads1.fastq.gz -pacbio-hifi reads2.fastq.gz -pacbio-hifi reads3.fastq.gz onSuccess="$base_path/on_init_complete.sh"
