#!/bin/bash
#Should be probably called within GraphAligner conda environment
#conda activate GraphAligner ; ./build.sh -j3 ; conda deactivate
set -eou

path=$(dirname $(readlink -e $0))

echo "Fetching submodules"
cd $path
git submodule update --init --recursive
cd -

echo "Building Canu"
cd $path/canu/src
make $@
cd -

echo "Building miniasm"
cd $path/miniasm
make $@
cd -

echo "Building gfacpp tool collection"
cd $path/gfacpp
make $@
cd -

echo "Building GraphAligner"
cd $path/GraphAligner
make $@
cd -
