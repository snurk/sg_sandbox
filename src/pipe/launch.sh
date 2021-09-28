#!/bin/bash
set -eou

if [ "$#" -lt 1 ]; then
    echo "script.sh <config.yaml> <out_dir> [additional snakemake arguments]*"
    exit
fi

dir=$(readlink -e $2)
config=$(readlink -e $1)

cd $(dirname $0)

export SCRIPT_PATH=$(readlink -e .)

# -o /home/snurk/logs/ -e /home/snurk/logs/
echo "Will use config file: \"$config\""
snakemake --latency-wait 60 -j 60 --local-cores 4 --cluster-config cluster.json --configfile $config --cluster "sbatch --ntasks 1 --mem {cluster.mem}G --cpus-per-task {cluster.n} --time {cluster.time}" --directory $dir "${@:3}"

cd -
