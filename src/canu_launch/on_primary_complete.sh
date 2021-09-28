#!/bin/bash
set -eou

if [ "$#" -lt 1 ]; then
    echo "script.sh <name>"
    exit 239
fi

name=$1
#absolute path
launch_path=$(pwd)

root_path=$launch_path/..
pipe_path=$root_path/microasm

#if [ ! -f $root_path/config.yaml ] ; then
#    echo "config.yaml file for running string-graph pipeline is missing from $root_path"
#fi

mkdir -p $pipe_path

sbatch --cpus-per-task=4 --mem=4g --time=24:00:00 --partition=norm --wrap "$(dirname $0)/../pipe/launch.sh $SG_CONFIG $pipe_path -C ASSEMBLY=\"$launch_path\""

#for debug
#sbatch --cpus-per-task=4 --mem=4g --time=2:00:00 --partition=norm --wrap "$(dirname $0)/../pipe/launch_local.sh $SG_CONFIG $pipe_path -C ASSEMBLY=\"$launch_path\""
#$(dirname $0)/../pipe/launch_local.sh $SG_CONFIG $pipe_path -C ASSEMBLY="$launch_path"
