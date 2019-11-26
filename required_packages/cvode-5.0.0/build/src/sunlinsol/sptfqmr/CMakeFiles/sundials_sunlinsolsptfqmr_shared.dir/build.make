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
include src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/depend.make

# Include the progress variables for this target.
include src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/progress.make

# Include the compile flags for this target's objects.
include src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/flags.make

src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/sunlinsol_sptfqmr.c.o: src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/flags.make
src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/sunlinsol_sptfqmr.c.o: ../src/sunlinsol/sptfqmr/sunlinsol_sptfqmr.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/sunlinsol_sptfqmr.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunlinsol/sptfqmr && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/sunlinsol_sptfqmr.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sunlinsol/sptfqmr/sunlinsol_sptfqmr.c

src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/sunlinsol_sptfqmr.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/sunlinsol_sptfqmr.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunlinsol/sptfqmr && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sunlinsol/sptfqmr/sunlinsol_sptfqmr.c > CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/sunlinsol_sptfqmr.c.i

src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/sunlinsol_sptfqmr.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/sunlinsol_sptfqmr.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunlinsol/sptfqmr && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sunlinsol/sptfqmr/sunlinsol_sptfqmr.c -o CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/sunlinsol_sptfqmr.c.s

src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_math.c.o: src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/flags.make
src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_math.c.o: ../src/sundials/sundials_math.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_math.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunlinsol/sptfqmr && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_math.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_math.c

src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_math.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_math.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunlinsol/sptfqmr && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_math.c > CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_math.c.i

src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_math.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_math.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunlinsol/sptfqmr && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_math.c -o CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_math.c.s

src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_nvector.c.o: src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/flags.make
src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_nvector.c.o: ../src/sundials/sundials_nvector.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_nvector.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunlinsol/sptfqmr && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_nvector.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nvector.c

src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_nvector.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_nvector.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunlinsol/sptfqmr && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nvector.c > CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_nvector.c.i

src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_nvector.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_nvector.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunlinsol/sptfqmr && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nvector.c -o CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_nvector.c.s

src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_linearsolver.c.o: src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/flags.make
src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_linearsolver.c.o: ../src/sundials/sundials_linearsolver.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_linearsolver.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunlinsol/sptfqmr && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_linearsolver.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_linearsolver.c

src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_linearsolver.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_linearsolver.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunlinsol/sptfqmr && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_linearsolver.c > CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_linearsolver.c.i

src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_linearsolver.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_linearsolver.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunlinsol/sptfqmr && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_linearsolver.c -o CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_linearsolver.c.s

src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_iterative.c.o: src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/flags.make
src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_iterative.c.o: ../src/sundials/sundials_iterative.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_iterative.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunlinsol/sptfqmr && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_iterative.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_iterative.c

src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_iterative.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_iterative.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunlinsol/sptfqmr && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_iterative.c > CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_iterative.c.i

src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_iterative.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_iterative.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunlinsol/sptfqmr && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_iterative.c -o CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_iterative.c.s

# Object files for target sundials_sunlinsolsptfqmr_shared
sundials_sunlinsolsptfqmr_shared_OBJECTS = \
"CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/sunlinsol_sptfqmr.c.o" \
"CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_math.c.o" \
"CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_nvector.c.o" \
"CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_linearsolver.c.o" \
"CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_iterative.c.o"

# External object files for target sundials_sunlinsolsptfqmr_shared
sundials_sunlinsolsptfqmr_shared_EXTERNAL_OBJECTS =

src/sunlinsol/sptfqmr/libsundials_sunlinsolsptfqmr.so.3.0.0: src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/sunlinsol_sptfqmr.c.o
src/sunlinsol/sptfqmr/libsundials_sunlinsolsptfqmr.so.3.0.0: src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_math.c.o
src/sunlinsol/sptfqmr/libsundials_sunlinsolsptfqmr.so.3.0.0: src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_nvector.c.o
src/sunlinsol/sptfqmr/libsundials_sunlinsolsptfqmr.so.3.0.0: src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_linearsolver.c.o
src/sunlinsol/sptfqmr/libsundials_sunlinsolsptfqmr.so.3.0.0: src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/__/__/sundials/sundials_iterative.c.o
src/sunlinsol/sptfqmr/libsundials_sunlinsolsptfqmr.so.3.0.0: src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/build.make
src/sunlinsol/sptfqmr/libsundials_sunlinsolsptfqmr.so.3.0.0: src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking C shared library libsundials_sunlinsolsptfqmr.so"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunlinsol/sptfqmr && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/link.txt --verbose=$(VERBOSE)
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunlinsol/sptfqmr && $(CMAKE_COMMAND) -E cmake_symlink_library libsundials_sunlinsolsptfqmr.so.3.0.0 libsundials_sunlinsolsptfqmr.so.3 libsundials_sunlinsolsptfqmr.so

src/sunlinsol/sptfqmr/libsundials_sunlinsolsptfqmr.so.3: src/sunlinsol/sptfqmr/libsundials_sunlinsolsptfqmr.so.3.0.0
	@$(CMAKE_COMMAND) -E touch_nocreate src/sunlinsol/sptfqmr/libsundials_sunlinsolsptfqmr.so.3

src/sunlinsol/sptfqmr/libsundials_sunlinsolsptfqmr.so: src/sunlinsol/sptfqmr/libsundials_sunlinsolsptfqmr.so.3.0.0
	@$(CMAKE_COMMAND) -E touch_nocreate src/sunlinsol/sptfqmr/libsundials_sunlinsolsptfqmr.so

# Rule to build all files generated by this target.
src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/build: src/sunlinsol/sptfqmr/libsundials_sunlinsolsptfqmr.so

.PHONY : src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/build

src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/clean:
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunlinsol/sptfqmr && $(CMAKE_COMMAND) -P CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/cmake_clean.cmake
.PHONY : src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/clean

src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/depend:
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0 /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sunlinsol/sptfqmr /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunlinsol/sptfqmr /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/sunlinsol/sptfqmr/CMakeFiles/sundials_sunlinsolsptfqmr_shared.dir/depend

