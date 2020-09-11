#!/bin/bash

if [ "$#" -lt 2 ]; then
    echo "script.sh <chr_name> <out_dir> [seed, default 239]"
    exit
fi

export CHR="$1"

dir=$(readlink -e $2)

if [ "$#" -lt 3 ]; then
    export RANDOM_SEED=239
else
    export RANDOM_SEED=$3
    echo "Random seed for layoutReads set at $RANDOM_SEED"
fi

cd $(dirname "$0")
export SCRIPT_PATH=$(readlink -e .)

echo "Launching consensus for $CHR"

# -o /home/snurk/logs/ -e /home/snurk/logs/
snakemake --latency-wait 60 -j 60 --local-cores 4 --cluster-config cluster.json --cluster "sbatch --ntasks 1 --mem {cluster.mem}G --cpus-per-task {cluster.n} --time {cluster.time}" --directory $dir &>> ${dir%/}.cns.log
# "${@:3}"

echo "Consensus for $CHR finished"
cd -
