#!/bin/bash
set -eou
#module load seqtk

if [ "$#" -lt 3 ]; then
    echo "script.sh <out_dir> <read trim len (non-hpc)> <reads>+"
    exit 239
fi

out_dir=$1
max_len=$2

mkdir -p $out_dir
rm -f $out_dir/reads.fa.gz
for reads in ${@:3} ; do
    seqtk trimfq -L $max_len $reads | seqtk seq -A | gzip --fast >> $out_dir/reads.fa.gz
done

base_path=$(dirname $(readlink -e $0))
export CANU_BIN=$(readlink -e $base_path/../../canu/build/bin/)
export MAX_HPC_LEN=1000000

$base_path/canu.sh $out_dir/canu $out_dir/reads.fa.gz onSuccess="$base_path/on_init_complete.sh"
