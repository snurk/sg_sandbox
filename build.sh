#!/bin/bash
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
