#!/usr/bin/env bash

#$HOME/local is the folder to which the libraries are being installed
if [ ! -d $HOME/local ]; then
    mkdir $HOME/local
fi

## this exports the location of the libraries
## It is very advisable to add these lines to the .bashrc file in your
## home directory, so that all these paths are available for each login at euler.
export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$HOME/local/lib:$LIBRARY_PATH
export CPATH=$HOME/local/include:$CPATH
export PATH=$HOME/local/bin:$PATH

######## installing muparser ######################
# This library is needed for the parsing of math formulas #
printf "\nInstalling muParser!\n"
cd ../required_packages/muparser-2.2.6.1/build # the location of the library
cmake ../ -DCMAKE_INSTALL_PREFIX=$HOME/local
make install
make clear
cd ../../

######## installing DPGMM #########################
# This library is needed for the density approximation#
printf "\nInstalling DP-GMM!\n"
cd required_packages/DP-GMM # the location of the library
if [ ! -d ./build ]; then
    mkdir ./build
fi
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/local/
make
make install
cd ../
rm -r build
cd ../../


######## installing cvode #########################
printf "\n\nInstalling cvode!\n"
cd required_packages/cvode-5.0.0
if [ ! -d ./build ]; then
    mkdir ./build
fi
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/local/
make
make install
cd ../
rm -r build
cd ../../