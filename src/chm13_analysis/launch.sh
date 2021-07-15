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

# -o /home/snurk/logs/ -e /home/snurk/logs/
snakemake --latency-wait 60 -j 60 --local-cores 4 --cluster-config cluster.json --cluster "sbatch --ntasks 1 --mem {cluster.mem}G --cpus-per-task {cluster.n} --time {cluster.time}" --directory $dir "${@:2}"

cd -
