#!/bin/bash
set -eou

path=$(dirname $(readlink -e $0))

echo "Cleaning Canu"
cd $path/canu/src
make clean
cd -

echo "Cleaning miniasm"
cd $path/miniasm
make clean
cd -

echo "Cleaning gfacpp tool collection"
cd $path/gfacpp
make clean
cd -
