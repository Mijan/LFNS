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
include examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/depend.make

# Include the progress variables for this target.
include examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/progress.make

# Include the compile flags for this target's objects.
include examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/flags.make

examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.o: examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/flags.make
examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.o: ../examples/cvode/serial/cvRoberts_dns_negsol.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.o"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/examples/cvode/serial && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.o   -c /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/examples/cvode/serial/cvRoberts_dns_negsol.c

examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.i"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/examples/cvode/serial && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/examples/cvode/serial/cvRoberts_dns_negsol.c > CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.i

examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.s"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/examples/cvode/serial && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/examples/cvode/serial/cvRoberts_dns_negsol.c -o CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.s

# Object files for target cvRoberts_dns_negsol
cvRoberts_dns_negsol_OBJECTS = \
"CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.o"

# External object files for target cvRoberts_dns_negsol
cvRoberts_dns_negsol_EXTERNAL_OBJECTS =

examples/cvode/serial/cvRoberts_dns_negsol: examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.o
examples/cvode/serial/cvRoberts_dns_negsol: examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/build.make
examples/cvode/serial/cvRoberts_dns_negsol: src/cvode/libsundials_cvode.so.5.0.0
examples/cvode/serial/cvRoberts_dns_negsol: src/nvector/serial/libsundials_nvecserial.so.5.0.0
examples/cvode/serial/cvRoberts_dns_negsol: /usr/lib/librt.so
examples/cvode/serial/cvRoberts_dns_negsol: examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable cvRoberts_dns_negsol"
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/examples/cvode/serial && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cvRoberts_dns_negsol.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/build: examples/cvode/serial/cvRoberts_dns_negsol

.PHONY : examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/build

examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/clean:
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/examples/cvode/serial && $(CMAKE_COMMAND) -P CMakeFiles/cvRoberts_dns_negsol.dir/cmake_clean.cmake
.PHONY : examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/clean

examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/depend:
	cd /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0 /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/examples/cvode/serial /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/examples/cvode/serial /home/jan/crypt/CLionProjects/LFNS/required_packages/cvode-5.0.0/build/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/depend

