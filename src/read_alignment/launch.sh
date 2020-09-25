#!/bin/bash

if [ "$#" -lt 1 ]; then
    echo "script.sh <out_dir> [additional snakemake arguments]*"
    exit
fi

dir=$(readlink -e $1)

cd $(dirname "$0")
export SCRIPT_PATH=$(readlink -e .)

export PATH=~/git/Winnowmap/bin/:~/git/Sniffles/bin/sniffles-core-1.0.12/:$PATH
module load samtools

# -o /home/snurk/logs/ -e /home/snurk/logs/
snakemake --latency-wait 60 -j 200 --local-cores 4 --cluster-config cluster.json --cluster "sbatch --ntasks 1 --mem {cluster.mem}G --cpus-per-task {cluster.n} --time {cluster.time}" --directory $dir "${@:2}"

cd -
