#!/bin/bash
set -eou

if [ "$#" -lt 1 ]; then
    echo "script.sh <name>"
    exit 239
fi

name=$1
launch_path=$(pwd)
oea_run_path=${launch_path}_oea

if [ ! -f oea_corrected.fa.gz ] ; then
    sbatch -W --cpus-per-task=8 --mem=20g --time=8:00:00 --partition=norm $(dirname $0)/fix_errors.sh ./ oea_corrected.fa.gz
fi

$(dirname $0)/canu.sh $oea_run_path ./oea_corrected.fa.gz
#onSuccess=$(dirname $0)/on_oea_run_complete.sh
