#!/usr/bin/env bash
## This script installs the simulation code into the subdirectory $HOME/LNSv4/local_euler_installs, which is assumed to exist allready.
## if it does not, make sure to run the 'install_libraries.sh' script first


## this load the boost library to euler
module load open_mpi boost/1.59.0 cmake/3.5.2
export CPPFLAGS="${CPPFLAGS} -I${BOOST_INCLUDEDIR}"
export LDFLAGS="-L${BOOST_LIBRARYDIR} ${LDFLAGS}"
export EIGEN3_ROOT=$HOME/LFNS/required_packages/eigen3.4

printf "\n\nInstalling Libraries\n"
source install_libraries.sh

printf "\n\nInstalling LFNS!\n"
source install_LFNS.sh


