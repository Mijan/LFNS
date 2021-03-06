# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2019, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# vector level CMakeLists.txt for SUNDIALS
# ------------------------------------------------------------------------------

# Always add SUNDIALS serial vector
add_subdirectory(serial)

# Always add SUNDIALS ManyVector
add_subdirectory(manyvector)

# Add NVECTOR modules with depedencies
if(MPI_ENABLE AND MPI_C_FOUND)
  add_subdirectory(parallel)
endif(MPI_ENABLE AND MPI_C_FOUND)

if(MPI_ENABLE AND MPI_C_FOUND)
  add_subdirectory(mpiplusx)
endif(MPI_ENABLE AND MPI_C_FOUND)

if(HYPRE_ENABLE AND HYPRE_FOUND)
  add_subdirectory(parhyp)
endif(HYPRE_ENABLE AND HYPRE_FOUND)

if(OPENMP_ENABLE AND OPENMP_FOUND)
  add_subdirectory(openmp)
endif(OPENMP_ENABLE AND OPENMP_FOUND)

if(OPENMP_DEVICE_ENABLE AND OPENMP_FOUND AND OPENMP_SUPPORTS_DEVICE_OFFLOADING)
  add_subdirectory(openmpdev)
endif(OPENMP_DEVICE_ENABLE AND OPENMP_FOUND AND OPENMP_SUPPORTS_DEVICE_OFFLOADING)

if(PTHREAD_ENABLE AND PTHREADS_FOUND)
  add_subdirectory(pthreads)
endif(PTHREAD_ENABLE AND PTHREADS_FOUND)

if(PETSC_ENABLE AND PETSC_FOUND)
  add_subdirectory(petsc)
endif(PETSC_ENABLE AND PETSC_FOUND)

if(CUDA_ENABLE AND CMAKE_CUDA_COMPILER)
  add_subdirectory(cuda)
endif(CUDA_ENABLE AND CMAKE_CUDA_COMPILER)

if(RAJA_ENABLE AND RAJA_FOUND)
  add_subdirectory(raja)
endif(RAJA_ENABLE AND RAJA_FOUND)

if(Trilinos_ENABLE AND Trilinos_FOUND)
  add_subdirectory(trilinos)
endif(Trilinos_ENABLE AND Trilinos_FOUND)
