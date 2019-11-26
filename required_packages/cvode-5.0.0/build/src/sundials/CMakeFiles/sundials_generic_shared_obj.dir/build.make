# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build

# Include any dependencies generated for this target.
include src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/depend.make

# Include the progress variables for this target.
include src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/progress.make

# Include the compile flags for this target's objects.
include src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/flags.make

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_band.c.o: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/flags.make
src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_band.c.o: ../src/sundials/sundials_band.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_band.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_band.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_band.c

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_band.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_generic_shared_obj.dir/sundials_band.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_band.c > CMakeFiles/sundials_generic_shared_obj.dir/sundials_band.c.i

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_band.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_generic_shared_obj.dir/sundials_band.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_band.c -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_band.c.s

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_dense.c.o: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/flags.make
src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_dense.c.o: ../src/sundials/sundials_dense.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_dense.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_dense.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_dense.c

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_dense.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_generic_shared_obj.dir/sundials_dense.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_dense.c > CMakeFiles/sundials_generic_shared_obj.dir/sundials_dense.c.i

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_dense.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_generic_shared_obj.dir/sundials_dense.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_dense.c -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_dense.c.s

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_direct.c.o: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/flags.make
src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_direct.c.o: ../src/sundials/sundials_direct.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_direct.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_direct.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_direct.c

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_direct.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_generic_shared_obj.dir/sundials_direct.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_direct.c > CMakeFiles/sundials_generic_shared_obj.dir/sundials_direct.c.i

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_direct.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_generic_shared_obj.dir/sundials_direct.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_direct.c -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_direct.c.s

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_iterative.c.o: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/flags.make
src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_iterative.c.o: ../src/sundials/sundials_iterative.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_iterative.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_iterative.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_iterative.c

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_iterative.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_generic_shared_obj.dir/sundials_iterative.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_iterative.c > CMakeFiles/sundials_generic_shared_obj.dir/sundials_iterative.c.i

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_iterative.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_generic_shared_obj.dir/sundials_iterative.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_iterative.c -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_iterative.c.s

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_linearsolver.c.o: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/flags.make
src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_linearsolver.c.o: ../src/sundials/sundials_linearsolver.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_linearsolver.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_linearsolver.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_linearsolver.c

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_linearsolver.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_generic_shared_obj.dir/sundials_linearsolver.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_linearsolver.c > CMakeFiles/sundials_generic_shared_obj.dir/sundials_linearsolver.c.i

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_linearsolver.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_generic_shared_obj.dir/sundials_linearsolver.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_linearsolver.c -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_linearsolver.c.s

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_math.c.o: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/flags.make
src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_math.c.o: ../src/sundials/sundials_math.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_math.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_math.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_math.c

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_math.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_generic_shared_obj.dir/sundials_math.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_math.c > CMakeFiles/sundials_generic_shared_obj.dir/sundials_math.c.i

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_math.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_generic_shared_obj.dir/sundials_math.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_math.c -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_math.c.s

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_matrix.c.o: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/flags.make
src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_matrix.c.o: ../src/sundials/sundials_matrix.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building C object src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_matrix.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_matrix.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_matrix.c

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_matrix.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_generic_shared_obj.dir/sundials_matrix.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_matrix.c > CMakeFiles/sundials_generic_shared_obj.dir/sundials_matrix.c.i

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_matrix.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_generic_shared_obj.dir/sundials_matrix.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_matrix.c -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_matrix.c.s

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_nonlinearsolver.c.o: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/flags.make
src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_nonlinearsolver.c.o: ../src/sundials/sundials_nonlinearsolver.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building C object src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_nonlinearsolver.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_nonlinearsolver.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nonlinearsolver.c

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_nonlinearsolver.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_generic_shared_obj.dir/sundials_nonlinearsolver.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nonlinearsolver.c > CMakeFiles/sundials_generic_shared_obj.dir/sundials_nonlinearsolver.c.i

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_nonlinearsolver.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_generic_shared_obj.dir/sundials_nonlinearsolver.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nonlinearsolver.c -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_nonlinearsolver.c.s

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector.c.o: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/flags.make
src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector.c.o: ../src/sundials/sundials_nvector.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building C object src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nvector.c

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nvector.c > CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector.c.i

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nvector.c -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector.c.s

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector_senswrapper.c.o: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/flags.make
src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector_senswrapper.c.o: ../src/sundials/sundials_nvector_senswrapper.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building C object src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector_senswrapper.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector_senswrapper.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nvector_senswrapper.c

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector_senswrapper.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector_senswrapper.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nvector_senswrapper.c > CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector_senswrapper.c.i

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector_senswrapper.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector_senswrapper.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nvector_senswrapper.c -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector_senswrapper.c.s

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_version.c.o: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/flags.make
src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_version.c.o: ../src/sundials/sundials_version.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building C object src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_version.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_version.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_version.c

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_version.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_generic_shared_obj.dir/sundials_version.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_version.c > CMakeFiles/sundials_generic_shared_obj.dir/sundials_version.c.i

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_version.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_generic_shared_obj.dir/sundials_version.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_version.c -o CMakeFiles/sundials_generic_shared_obj.dir/sundials_version.c.s

sundials_generic_shared_obj: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_band.c.o
sundials_generic_shared_obj: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_dense.c.o
sundials_generic_shared_obj: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_direct.c.o
sundials_generic_shared_obj: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_iterative.c.o
sundials_generic_shared_obj: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_linearsolver.c.o
sundials_generic_shared_obj: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_math.c.o
sundials_generic_shared_obj: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_matrix.c.o
sundials_generic_shared_obj: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_nonlinearsolver.c.o
sundials_generic_shared_obj: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector.c.o
sundials_generic_shared_obj: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_nvector_senswrapper.c.o
sundials_generic_shared_obj: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/sundials_version.c.o
sundials_generic_shared_obj: src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/build.make

.PHONY : sundials_generic_shared_obj

# Rule to build all files generated by this target.
src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/build: sundials_generic_shared_obj

.PHONY : src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/build

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/clean:
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials && $(CMAKE_COMMAND) -P CMakeFiles/sundials_generic_shared_obj.dir/cmake_clean.cmake
.PHONY : src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/clean

src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/depend:
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0 /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/sundials/CMakeFiles/sundials_generic_shared_obj.dir/depend

