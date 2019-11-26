# -----------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
#                David J. Gardner @ LLNL
# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------
# CMakeLists.txt for CVODE serial examples
#
# This file is generated from a template using variables
# set at configuration time. It can be used as a template for
# other user CMakeLists configuration files.
# -----------------------------------------------------------------

cmake_minimum_required(VERSION 3.5)

# Specify project name
project(CVODE_serial_examples C)

# Enable testing
include(CTest)

mark_as_advanced(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)

# Examples using SUNDIALS linear solvers
set(examples  cvRoberts_dns cvRoberts_dns_uw cvRoberts_dns_negsol cvRoberts_dns_constraints cvAdvDiff_bnd cvDirectDemo_ls cvDiurnal_kry_bp cvDiurnal_kry cvDisc_dns cvKrylovDemo_ls cvKrylovDemo_prec)
set(examples_dependencies )

# Examples using LAPACK linear solvers
set(examples_bl )
set(examples_dependencies_bl )

# Examples using KLU linear solvers
set(examples_klu )
set(examples_dependencies_klu )

# Examples using SuperLU_MT linear solvers
set(examples_slumt )
set(examples_dependencies_slumt )

# Specify path to SUNDIALS header files
set(SUNDIALS_INC_DIR
  /usr/local/include
  CACHE STRING
  "Location of SUNDIALS header files")

# Add path to SUNDIALS header files
include_directories(${SUNDIALS_INC_DIR})

# Set search path for SUNDIALS libraries
set(SUNDIALS_LIB_DIR /usr/local/lib)

# Find the SUNDIALS solver's library
find_library(SUNDIALS_SOLVER_LIB
  sundials_cvode ${SUNDIALS_LIB_DIR}
  DOC "CVODE library")

# Find the NVECTOR library
find_library(SUNDIALS_NVEC_LIB
  sundials_nvecserial ${SUNDIALS_LIB_DIR}
  DOC "NVECTOR library")

# Find the NVECTOR ManyVector library
find_library(SUNDIALS_MANYVEC_LIB
  sundials_nvecmanyvector ${SUNDIALS_LIB_DIR}
  DOC "NVECTOR manyvector library")

# Set an extra link directory if necessary (false if empty)
if(EXISTS )
  link_directories()
endif()

# Set additional libraries
set(SUNDIALS_EXTRA_LIB  -lm /usr/lib/librt.so CACHE STRING "Additional libraries")

# Set other libraries libraries
set(LAPACK_LIBRARIES  CACHE STRING "Lapack libraries")

include_directories()
set(SUPERLUMT_LIBRARIES   CACHE STRING "SuperLUMT libraries")
set(SUPERLUMT_THREAD_TYPE  CACHE STRING "SuperLUMT Threading")

include_directories()
set(KLU_LIBRARIES  CACHE STRING "KLU libraries")

# List of Sundials libraries shared across all examples
set(SUNDIALS_LIBS ${SUNDIALS_SOLVER_LIB} ${SUNDIALS_NVEC_LIB} ${SUNDIALS_MANYVEC_LIB} ${SUNDIALS_EXTRA_LIB})

# Build examples with SUNDIALS linear solvers one by one
foreach(example ${examples})

  # example source files
  add_executable(${example} ${example}.c ${examples_dependencies})

  # libraries to link against
  target_link_libraries(${example} ${SUNDIALS_LIBS})

  # add the example to ctest
  add_test(NAME ${example} COMMAND ${example})

endforeach()

# Build each example with LAPACK linear solvers one by one
if(LAPACK_LIBRARIES)
  foreach(example ${examples_bl})

    # Find the Sundials linear solver libraries
    find_library(SUNDIALS_SUNLINSOLLAPACKBAND_LIB
      sundials_sunlinsollapackband ${SUNDIALS_LIB_DIR}
      DOC "Sundials banded LAPACK linear solver library")

    find_library(SUNDIALS_SUNLINSOLLAPACKDENSE_LIB
      sundials_sunlinsollapackdense ${SUNDIALS_LIB_DIR}
      DOC "Sundials dense LAPACK linear solver library")

    # example source files
    add_executable(${example} ${example}.c ${examples_dependencies_bl})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS})
    target_link_libraries(${example} ${SUNDIALS_SUNLINSOLLAPACKBAND_LIB})
    target_link_libraries(${example} ${SUNDIALS_SUNLINSOLLAPACKDENSE_LIB})
    target_link_libraries(${example} ${LAPACK_LIBRARIES})

    # add the example to ctest
    add_test(NAME ${example} COMMAND ${example})

  endforeach(example ${examples_bl})
endif()


# Build each example with KLU linear solvers one by one
if(KLU_LIBRARIES)
  foreach(example ${examples_klu})

    # Find the Sundials linear solver libraries
    find_library(SUNDIALS_SUNLINSOLKLU_LIB
      sundials_sunlinsolklu ${SUNDIALS_LIB_DIR}
      DOC "Sundials KLU linear solver library")

    # example source files
    add_executable(${example} ${example}.c ${examples_dependencies_klu})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS})
    target_link_libraries(${example} ${SUNDIALS_SUNLINSOLKLU_LIB})
    target_link_libraries(${example} ${KLU_LIBRARIES})

    # add the example to ctest
    add_test(NAME ${example} COMMAND ${example})

  endforeach(example ${examples_klu})
endif()


# Build each example with SuperLU_MT linear solvers one by one
if(SUPERLUMT_LIBRARIES)

  # Find the Sundials linear solver libraries
  find_library(SUNDIALS_SUNLINSOLSLUMT_LIB
    sundials_sunlinsolsuperlumt ${SUNDIALS_LIB_DIR}
    DOC "Sundials KLU linear solver library")

  if(SUPERLUMT_THREAD_TYPE STREQUAL "OPENMP")
    find_package(OpenMP REQUIRED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  else(SUPERLUMT_THREAD_TYPE STREQUAL "PTHREAD")
    find_package(Threads REQUIRED)
  endif()

  foreach(example ${examples_slumt})

    # example source files
    add_executable(${example} ${example}.c ${examples_dependencies_slumt})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS})
    target_link_libraries(${example} ${SUNDIALS_SUNLINSOLSLUMT_LIB})
    target_link_libraries(${example} ${SUPERLUMT_LIBRARIES})

    if(SUPERLUMT_THREAD_TYPE STREQUAL "PTHREAD")
      target_link_libraries(${example} Threads::Threads)
    endif()

    # add the example to ctest
    add_test(NAME ${example} COMMAND ${example})

  endforeach(example ${examples_slumt})
endif()
