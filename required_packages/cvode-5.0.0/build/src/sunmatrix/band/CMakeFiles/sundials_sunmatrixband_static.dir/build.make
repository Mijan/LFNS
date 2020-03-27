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
include src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/depend.make

# Include the progress variables for this target.
include src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/progress.make

# Include the compile flags for this target's objects.
include src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/flags.make

src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/sunmatrix_band.c.o: src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/flags.make
src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/sunmatrix_band.c.o: ../src/sunmatrix/band/sunmatrix_band.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/sunmatrix_band.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunmatrix/band && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_sunmatrixband_static.dir/sunmatrix_band.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sunmatrix/band/sunmatrix_band.c

src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/sunmatrix_band.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_sunmatrixband_static.dir/sunmatrix_band.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunmatrix/band && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sunmatrix/band/sunmatrix_band.c > CMakeFiles/sundials_sunmatrixband_static.dir/sunmatrix_band.c.i

src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/sunmatrix_band.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_sunmatrixband_static.dir/sunmatrix_band.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunmatrix/band && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sunmatrix/band/sunmatrix_band.c -o CMakeFiles/sundials_sunmatrixband_static.dir/sunmatrix_band.c.s

src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_nvector.c.o: src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/flags.make
src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_nvector.c.o: ../src/sundials/sundials_nvector.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_nvector.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunmatrix/band && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_nvector.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nvector.c

src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_nvector.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_nvector.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunmatrix/band && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nvector.c > CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_nvector.c.i

src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_nvector.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_nvector.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunmatrix/band && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nvector.c -o CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_nvector.c.s

src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_matrix.c.o: src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/flags.make
src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_matrix.c.o: ../src/sundials/sundials_matrix.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_matrix.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunmatrix/band && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_matrix.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_matrix.c

src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_matrix.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_matrix.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunmatrix/band && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_matrix.c > CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_matrix.c.i

src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_matrix.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_matrix.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunmatrix/band && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_matrix.c -o CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_matrix.c.s

src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_math.c.o: src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/flags.make
src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_math.c.o: ../src/sundials/sundials_math.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_math.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunmatrix/band && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_math.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_math.c

src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_math.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_math.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunmatrix/band && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_math.c > CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_math.c.i

src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_math.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_math.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunmatrix/band && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_math.c -o CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_math.c.s

# Object files for target sundials_sunmatrixband_static
sundials_sunmatrixband_static_OBJECTS = \
"CMakeFiles/sundials_sunmatrixband_static.dir/sunmatrix_band.c.o" \
"CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_nvector.c.o" \
"CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_matrix.c.o" \
"CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_math.c.o"

# External object files for target sundials_sunmatrixband_static
sundials_sunmatrixband_static_EXTERNAL_OBJECTS =

src/sunmatrix/band/libsundials_sunmatrixband.a: src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/sunmatrix_band.c.o
src/sunmatrix/band/libsundials_sunmatrixband.a: src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_nvector.c.o
src/sunmatrix/band/libsundials_sunmatrixband.a: src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_matrix.c.o
src/sunmatrix/band/libsundials_sunmatrixband.a: src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/__/__/sundials/sundials_math.c.o
src/sunmatrix/band/libsundials_sunmatrixband.a: src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/build.make
src/sunmatrix/band/libsundials_sunmatrixband.a: src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking C static library libsundials_sunmatrixband.a"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunmatrix/band && $(CMAKE_COMMAND) -P CMakeFiles/sundials_sunmatrixband_static.dir/cmake_clean_target.cmake
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunmatrix/band && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sundials_sunmatrixband_static.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/build: src/sunmatrix/band/libsundials_sunmatrixband.a

.PHONY : src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/build

src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/clean:
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunmatrix/band && $(CMAKE_COMMAND) -P CMakeFiles/sundials_sunmatrixband_static.dir/cmake_clean.cmake
.PHONY : src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/clean

src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/depend:
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0 /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sunmatrix/band /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunmatrix/band /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/sunmatrix/band/CMakeFiles/sundials_sunmatrixband_static.dir/depend
