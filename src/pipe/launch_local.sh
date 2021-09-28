#!/bin/bash
set -eou

if [ "$#" -lt 1 ]; then
    echo "script.sh <with config.yaml> <out_dir> [additional snakemake arguments]*"
    exit
fi

config=$(readlink -e $1)
dir=$(readlink -e $2)

cd $(dirname $0)

export SCRIPT_PATH=$(readlink -e .)

snakemake --latency-wait 60 -j 8 --configfile $config --directory $dir "${@:3}"

cd -
