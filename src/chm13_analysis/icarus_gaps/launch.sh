#!/bin/bash

if [ "$#" -lt 1 ]; then
    echo "script.sh <out_dir> [additional snakemake arguments]*"
    exit
fi

module load quast
module load seqtk

snakemake --latency-wait 60 -j 60 --local-cores 4 --cluster "sbatch --ntasks 1 --mem 8g --cpus-per-task 8 --time 4:00:00" --directory "$@"
