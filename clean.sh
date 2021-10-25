#!/bin/bash
set -eou

path=$(dirname $(readlink -e $0))

echo "Building Canu"
cd $path/canu/src
make clean
cd -

echo "Building miniasm"
cd $path/miniasm
make clean
cd -

echo "Building gfacpp tool collection"
cd $path/gfacpp
make clean
cd -
