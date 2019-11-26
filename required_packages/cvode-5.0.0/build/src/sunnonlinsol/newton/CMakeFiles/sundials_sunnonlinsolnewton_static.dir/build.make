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
include src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/depend.make

# Include the progress variables for this target.
include src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/progress.make

# Include the compile flags for this target's objects.
include src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/flags.make

src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/sunnonlinsol_newton.c.o: src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/flags.make
src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/sunnonlinsol_newton.c.o: ../src/sunnonlinsol/newton/sunnonlinsol_newton.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/sunnonlinsol_newton.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunnonlinsol/newton && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_sunnonlinsolnewton_static.dir/sunnonlinsol_newton.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sunnonlinsol/newton/sunnonlinsol_newton.c

src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/sunnonlinsol_newton.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_sunnonlinsolnewton_static.dir/sunnonlinsol_newton.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunnonlinsol/newton && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sunnonlinsol/newton/sunnonlinsol_newton.c > CMakeFiles/sundials_sunnonlinsolnewton_static.dir/sunnonlinsol_newton.c.i

src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/sunnonlinsol_newton.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_sunnonlinsolnewton_static.dir/sunnonlinsol_newton.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunnonlinsol/newton && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sunnonlinsol/newton/sunnonlinsol_newton.c -o CMakeFiles/sundials_sunnonlinsolnewton_static.dir/sunnonlinsol_newton.c.s

src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_math.c.o: src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/flags.make
src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_math.c.o: ../src/sundials/sundials_math.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_math.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunnonlinsol/newton && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_math.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_math.c

src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_math.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_math.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunnonlinsol/newton && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_math.c > CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_math.c.i

src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_math.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_math.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunnonlinsol/newton && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_math.c -o CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_math.c.s

src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector.c.o: src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/flags.make
src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector.c.o: ../src/sundials/sundials_nvector.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunnonlinsol/newton && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nvector.c

src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunnonlinsol/newton && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nvector.c > CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector.c.i

src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunnonlinsol/newton && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nvector.c -o CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector.c.s

src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector_senswrapper.c.o: src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/flags.make
src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector_senswrapper.c.o: ../src/sundials/sundials_nvector_senswrapper.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector_senswrapper.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunnonlinsol/newton && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector_senswrapper.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nvector_senswrapper.c

src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector_senswrapper.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector_senswrapper.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunnonlinsol/newton && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nvector_senswrapper.c > CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector_senswrapper.c.i

src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector_senswrapper.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector_senswrapper.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunnonlinsol/newton && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nvector_senswrapper.c -o CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector_senswrapper.c.s

src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nonlinearsolver.c.o: src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/flags.make
src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nonlinearsolver.c.o: ../src/sundials/sundials_nonlinearsolver.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nonlinearsolver.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunnonlinsol/newton && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nonlinearsolver.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nonlinearsolver.c

src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nonlinearsolver.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nonlinearsolver.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunnonlinsol/newton && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nonlinearsolver.c > CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nonlinearsolver.c.i

src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nonlinearsolver.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nonlinearsolver.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunnonlinsol/newton && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sundials/sundials_nonlinearsolver.c -o CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nonlinearsolver.c.s

# Object files for target sundials_sunnonlinsolnewton_static
sundials_sunnonlinsolnewton_static_OBJECTS = \
"CMakeFiles/sundials_sunnonlinsolnewton_static.dir/sunnonlinsol_newton.c.o" \
"CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_math.c.o" \
"CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector.c.o" \
"CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector_senswrapper.c.o" \
"CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nonlinearsolver.c.o"

# External object files for target sundials_sunnonlinsolnewton_static
sundials_sunnonlinsolnewton_static_EXTERNAL_OBJECTS =

src/sunnonlinsol/newton/libsundials_sunnonlinsolnewton.a: src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/sunnonlinsol_newton.c.o
src/sunnonlinsol/newton/libsundials_sunnonlinsolnewton.a: src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_math.c.o
src/sunnonlinsol/newton/libsundials_sunnonlinsolnewton.a: src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector.c.o
src/sunnonlinsol/newton/libsundials_sunnonlinsolnewton.a: src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nvector_senswrapper.c.o
src/sunnonlinsol/newton/libsundials_sunnonlinsolnewton.a: src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/__/__/sundials/sundials_nonlinearsolver.c.o
src/sunnonlinsol/newton/libsundials_sunnonlinsolnewton.a: src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/build.make
src/sunnonlinsol/newton/libsundials_sunnonlinsolnewton.a: src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking C static library libsundials_sunnonlinsolnewton.a"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunnonlinsol/newton && $(CMAKE_COMMAND) -P CMakeFiles/sundials_sunnonlinsolnewton_static.dir/cmake_clean_target.cmake
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunnonlinsol/newton && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sundials_sunnonlinsolnewton_static.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/build: src/sunnonlinsol/newton/libsundials_sunnonlinsolnewton.a

.PHONY : src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/build

src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/clean:
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunnonlinsol/newton && $(CMAKE_COMMAND) -P CMakeFiles/sundials_sunnonlinsolnewton_static.dir/cmake_clean.cmake
.PHONY : src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/clean

src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/depend:
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0 /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/src/sunnonlinsol/newton /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunnonlinsol/newton /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/sunnonlinsol/newton/CMakeFiles/sundials_sunnonlinsolnewton_static.dir/depend

