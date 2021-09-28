#!/bin/bash
set -eou
#module load seqtk

if [ "$#" -lt 4 ]; then
    echo "script.sh <SG pipeline config (YAML)> <out_dir> <read trim len (non-hpc)> <reads>+"
    echo "Don't forget to satisfy all SG pipeline dependancies in the running environment!"
    exit 239
fi

config=$1
out_dir=$2
max_len=$3

mkdir -p $out_dir

if [ ! -f $out_dir/reads.fa.gz ] ; then
    rm -f $out_dir/reads.fa.gz.tmp
    for reads in ${@:4} ; do
        seqtk trimfq -L $max_len $reads | seqtk seq -A | gzip --fast >> $out_dir/reads.fa.gz.tmp
    done
    mv $out_dir/reads.fa.gz.tmp $out_dir/reads.fa.gz
fi

base_path=$(dirname $(readlink -e $0))
export CANU_BIN=$(readlink -e $base_path/../../canu/build/bin/)
#Threshold for ignoring reads in Bogart output
export MAX_BOGART_LEN=1000000
export SG_CONFIG=$(readlink -e $config)

$base_path/canu.sh $out_dir/canu $out_dir/reads.fa.gz onSuccess="$base_path/on_init_complete.sh"
