#!/bin/bash
set -eou

#OBSOLETE!!! Use master_trim.sh instead!

base_path=$(dirname $(readlink -e $0))
export CANU_BIN=$(readlink -e $base_path/../../canu/build/bin/)
export MAX_BOGART_LEN=18000

#example 1 -- single read file
$base_path/canu.sh out_canu reads.fa.gz onSuccess="$base_path/on_init_complete.sh"

#example 2 -- multiple read files
$base_path/canu.sh out_canu reads1.fastq.gz -pacbio-hifi reads2.fastq.gz -pacbio-hifi reads3.fastq.gz onSuccess="$base_path/on_init_complete.sh"
