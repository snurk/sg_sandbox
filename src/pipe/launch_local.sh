#!/bin/bash
set -eou

if [ "$#" -lt 1 ]; then
    echo "script.sh <out_dir (with config.yaml)> [additional snakemake arguments]*"
    exit
fi

module load snakemake
dir=$(readlink -e $1)

cd $(dirname $0)

export SCRIPT_PATH=$(readlink -e .)

snakemake --latency-wait 60 -j 8 --directory $dir "${@:2}"

cd -
