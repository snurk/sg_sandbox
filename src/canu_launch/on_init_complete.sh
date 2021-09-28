#!/bin/bash
set -eou

if [ "$#" -lt 1 ]; then
    echo "script.sh <name>"
    exit 239
fi

base_path=$(dirname $(readlink -e $0))

name=$1
launch_path=$(pwd)
#Absolute path
oea_run_path=${launch_path}_oea

if [ ! -f oea_corrected.fa.gz ] ; then
    sbatch -W --cpus-per-task=8 --mem=20g --time=8:00:00 --partition=norm $(dirname $0)/fix_errors.sh ./ oea_corrected.fa.gz
fi

$(dirname $0)/canu.sh $oea_run_path ./oea_corrected.fa.gz onSuccess="$base_path/on_primary_complete.sh"
