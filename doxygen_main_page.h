//
// Created by jan on 3/20/20.
//


/*! \mainpage The LFNS Toolbox
 *
 * \section General Information
 * \author Jan Mikelson <jan.mikelson@bsse.ethz.ch>
 *
 * This toolbox provides mainly and implementation of the likelihood-free Nested Sampling algorithm as described in
 * Mikelson, Jan, and Mustafa Khammash. "Likelihood-free nested sampling for biochemical reaction networks." bioRxiv (2019): 564047.
 *
 * During the development of the LFNS algorithm several other features had to be implemented as well, such as the simulation of deterministic,
 * stochastic and hybrid models of chemical reaction networks, particle filters and various sampling schemes. All these subroutines and methods are
 * collected in this toolbox.
 *
 * \section Installation
 * \subsection Required Libraries
 *
 * The LFNS toolbox depends on several packages that are all provided in the <a href="https://github.com/Mijan/LFNS/tree/publishable/required_packages">required_packages</a> folder,
 * but can also be downloaded from their corresponding websites
 *  - The <a href="http://eigen.tuxfamily.org/index.php?title=Main_Page">Eigen</a> library. Very powerful linear algebra library that provides very convenient matrix and vector classes
 *  - The <a href="https://beltoforion.de/article.php?a=muparser">muParser</a> library. A parsing library that simplifies parsing mathematical formulas from text files.
 *  - The <a href="https://www.boost.org/">boost</a> library. General utility library.
 *  - The <a href="https://github.com/Mijan/DPGMM">DP-GMM</a> library. A library that implements density estimation through a Gaussian Mixture.
 *  - The <a href="https://computing.llnl.gov/projects/sundials/cvode">cvode</a> library.
 *
 * \subsection Running cmake
 *
 * The installation is based on cmake, which means cmake needs to be installed <a href="https://cmake.org/'>https://cmake.org/</a>.
 * To install the LFNS-toolbox, a c++ compiler is required and the above mentioned libraries need to be installed and linked. In the following we give an example of
 * an installation procedure on a linux system.
 *
 * First, we pick a folder to which the toolbox and the required libraries are to be installed. In this example we pick the
 * location to be /usr/local. Note that on most linux systems this folder will already exists and will be appropriately linked,
 * so the following steps will most likely not be necessary.
 *  - First create the folder, including subfolders for libraries, header files and executables. Some of these commands may require super user privileges (use these commands with "sudo")
 *   \code{.sh}
 *  sudo mkdir /usr/local /usr/local/lib /usr/local/include /usr/local/bin
 *  \endcode
 *
 *  - Make sure this folder is linked whenever a new termina is opemend. For this open the .bashrc file in the home directory and add the following lines.
 *  \code{.sh}
 *  export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
 *  export LIBRARY_PATH=/usr/local/lib:$LIBRARY_PATH
 *  export CPATH=/usr/local/include:$CPATH
 *  export PATH=/usr/local/bin:$PATH
 *  \endcode
 *  -Then install all required libraries. For this assume that the folder containing the LFNS toolbox is in ~/LFNS. The libraries are all installed using cmake.
 *  For this change to each of the build subdirectories and call cmake and make install
 * \code{.sh}
 *  cd ~/LFNS/requred_packages/cvode-5.0.0/build
 *  cmake ../ -DCMAKE_INSTALL_PREFIX=/usr/local/
 *  make install
 *  \endcode
 * \code{.sh}
 *  cd ~/LFNS/requred_packages/DP-GMM/build
 *  cmake ../ -DCMAKE_INSTALL_PREFIX=/usr/local/
 *  make install
 *  \endcode
 * \code{.sh}
 *  cd ~/LFNS/requred_packages/muparser-2.2.6.1/build
 *  cmake ../ -DCMAKE_INSTALL_PREFIX=/usr/local/
 *  make install
 *  \endcode
 *  Eigen is a pure header library and can be either installed through the systems package manager (for Ubuntu for instance
 *  "apt-get install eigen") or the content of the eigen.34 folder can just be copied in the /usr/local/include directory.
 *  - After all these libraries have been installed, the LFNS toolbox can be installed by switching into the build directory and calling again cmake and make install
 * \code{.sh}
 *  cd ~/LFNS/build
 *  cmake ../ -DCMAKE_INSTALL_PREFIX=/usr/local/
 *  make install
 *  \endcode
 *  To check if everything has been installed correctly, try typing \code{.sh} lfns_seq --help \endcode
 *  This should produce a help message.
 */
