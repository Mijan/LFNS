//
// Created by jan on 3/20/20.
//


/*! \mainpage The LFNS Toolbox
 *
 * \section General Information
 * \author Jan Mikelson <jan.mikelson@bsse.ethz.ch>
 *
 * This toolbox provides mainly an implementation of the likelihood-free Nested Sampling algorithm as described in
 * Mikelson, Jan, and Mustafa Khammash. "Likelihood-free nested sampling for biochemical reaction networks." bioRxiv (2019): 564047.
 *
 * During the development of the LFNS algorithm several other features had to be implemented as well, such as the simulation of deterministic,
 * stochastic and hybrid models of chemical reaction networks, particle filters and various sampling schemes. All these subroutines and methods are
 * collected in this toolbox.
 *
 * \section Installation
 * \subsection Required_Libraries
 *
 * The LFNS toolbox depends on several packages that are all provided in the <a href="https://github.com/Mijan/LFNS/tree/publishable/required_packages">required_packages</a> folder,
 * but can also be downloaded from their corresponding websites
 *  - The <a href="http://eigen.tuxfamily.org/index.php?title=Main_Page">Eigen</a> library. Very powerful linear algebra library that provides very convenient matrix and vector classes
 *  - The <a href="https://beltoforion.de/article.php?a=muparser">muParser</a> library. A parsing library that simplifies parsing mathematical formulas from text files.
 *  - The <a href="https://www.boost.org/">boost</a> library. General utility library.
 *  - The <a href="https://github.com/Mijan/DPGMM">DP-GMM</a> library. A library that implements density estimation through a Gaussian Mixture.
 *  - The <a href="https://computing.llnl.gov/projects/sundials/cvode">cvode</a> library.
 *
 * \subsection setting_up_folder Setting up installation folder
 *
 * The installation is based on cmake, which means cmake needs to be installed <a href="https://cmake.org/">https://cmake.org/</a>.
 * To install the LFNS-toolbox, a c++ compiler is required and the above mentioned libraries need to be installed and linked. In the following we give an example of
 * an installation procedure on a linux system.
 *
 * First, we pick a folder to which the toolbox and the required libraries are to be installed. We will denote this folder with /home/install_dest. On Linux based systems we would usually use the folder /usr/local.
 * On most Linux distributions this folder is already existent and properly linked. I
 *  - In the case where the install folder does not exist yet, we create the folder including subfolders for libraries, header files and executables. So assuming we want to install the toolbox into the folder /home/install_dest we type
 *   \code{.sh}
 *  mkdir  /home/install_dest/lib /home/install_dest/include /home/install_dest/bin
 *  \endcode
 *
 *  - Make sure this folder is linked whenever a new terminal is opened. For this open the .bashrc file in the home directory and add the following lines.
 *  \code{.sh}
 *  export LD_LIBRARY_PATH=/home/install_dest/lib:$LD_LIBRARY_PATH
 *  export LIBRARY_PATH=/home/install_dest/lib:$LIBRARY_PATH
 *  export CPATH=/home/install_dest/include:$CPATH
 *  export PATH=/home/install_dest/bin:$PATH
 *  \endcode
 *
 * \subsection install_required_libs Installing required libraries
 *  - The Boost library is a standard utility library and can be installed using the standard package manager (for Ubuntu it would be "apt-get install boost")
 *  or downloaded and install from the official website <a href="https://www.boost.org/doc/libs/1_61_0/more/getting_started/unix-variants.html">here</a>
 *  - Then install all required libraries. For this assume that the folder containing the LFNS toolbox is in /home/LFNS. The libraries are all installed using cmake.
 *  For this change to each of the build subdirectories and call cmake and make install
 * \code{.sh}
 *  cd /home/LFNS/requred_packages/cvode-5.0.0/build
 *  cmake ../ -DCMAKE_INSTALL_PREFIX=/home/install_dest/
 *  make install
 *  \endcode
 * \code{.sh}
 *  cd /home/LFNS/requred_packages/DP-GMM/build
 *  cmake ../ -DCMAKE_INSTALL_PREFIX=/home/install_dest/
 *  make install
 *  \endcode
 * \code{.sh}
 *  cd /home/LFNS/requred_packages/muparser-2.2.6.1/build
 *  cmake ../ -DCMAKE_INSTALL_PREFIX=/home/install_dest/
 *  make install
 *  \endcode
 *  Eigen is a pure header library and can be either installed through the systems package manager (for Ubuntu for instance
 *  "apt-get install eigen") or the content of the eigen.34 folder can just be copied in the /usr/local/include directory.
 *
 *  \subsection  install_lfns Installing the LFNS toolbox
 *  - After all these libraries have been installed, the LFNS toolbox can be installed by switching into the build directory and calling again cmake and make install
 * \code{.sh}
 *  cd /home/LFNS/build
 *  cmake ../ -DCMAKE_INSTALL_PREFIX=/home/install_dest
 *  make install
 *  \endcode
 *  To check if everything has been installed correctly, try typing \code{.sh} lfns_seq --help \endcode
 *  This should produce a help message.
 *
 *  \section running_lfns Using the toolbox
 *  Currently the toolbox provides four callable functions:
 *      - \a lfns_seq runs the LFNS algorithm sequentially (i.e. using only one thread)
 *      - \a lfns_mpi runs the LFNS algorithm in parallel, using MPI.
 *      - \a simulation runs a simulation algorithm to simulate a ODE system, run the SSA algorithm or perform a hybrid simulation
 *      - \a compute_likelihood runs a particle filter to estimate the likelihood of a particular observation for a given model and a parameter vector
 *
 *  In the following we will briefly outline how to use these commands. All these commands have in common that they take a configuration file as the first arguments, containing all the necessary information about the model under ocnsideration.
 *
 *  \subsection config_file Configuration file
 *
 *  To run any of the above mentioned commands we need to provide a configuration file containing all the required information about the model, the data, parameter priors and so on.
 *  We will not describe the configuration file here in all detail, but will rather briefly outline their structure and refer to the examples <a href="https://github.com/Mijan/LFNS/tree/publishable/examples">here</a>.
 *
 *
 */
