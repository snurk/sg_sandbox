#!/bin/bash
set -eou

if [ "$#" -lt 3 ]; then
    echo "script.sh <chr_name> <config.yaml> <out_dir> [seed, default 239]"
    exit 239
fi

export CHR="$1"
config=$2
dir=$3

if [ "$#" -lt 4 ]; then
    export RANDOM_SEED=239
else
    export RANDOM_SEED=$4
    echo "Random seed for layoutReads set at $RANDOM_SEED"
fi

if [ ! -f $config ] ; then
    echo "Config file $config doesn't exist"
    exit 239
fi

if [ ! -f $dir/layout.txt ] ; then
    echo "layout.txt file wasn't present in the output folder ($dir)"
    exit 239
fi

config=$(readlink -e $config)
dir=$(readlink -e $dir)

export CANU_BIN=$(dirname $(readlink -e $0))/../../canu/build/bin/

cd $(dirname "$0")
export SCRIPT_PATH=$(readlink -e .)

echo "Launching consensus for $CHR; output in $dir"

# -o /home/snurk/logs/ -e /home/snurk/logs/
snakemake --latency-wait 60 -j 60 --local-cores 4 --cluster-config cluster.json --configfile $config --cluster "sbatch --ntasks 1 --mem {cluster.mem}G --cpus-per-task {cluster.n} --time {cluster.time}" --directory $dir #&>> ${dir%/}.cns.log
# "${@:3}"

echo "Consensus for $CHR finished"
cd -
