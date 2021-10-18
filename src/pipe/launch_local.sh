#!/bin/bash
set -eou

if [ "$#" -lt 2 ]; then
    echo "script.sh <config.yaml> <out_dir> [additional snakemake arguments]*"
    exit 239
fi

config=$(readlink -e $1)
mkdir -p $2
dir=$(readlink -e $2)

cd $(dirname $0)

export SCRIPT_PATH=$(readlink -e .)

echo "Will use config file: \"$config\""
echo "Output folder: \"$dir\""
echo "Additional snakemake arguments: \"${@:3}\""

snakemake --latency-wait 60 -j 8 --configfile $config --directory $dir "${@:3}"

cd -
