#!/usr/bin/env bash

cd ../
## installing the toolbox
if [ ! -d ./build ]; then
    mkdir ./build
fi
cd build
printf "\nCalling CMake\n"
export EIGEN3_ROOT=$HOME/LFNS/required_packages/eigen3.4
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/local/
printf "\nCalling Make\n"
make
make install
cd ../