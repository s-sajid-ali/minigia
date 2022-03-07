# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/sajid/packages/spack/opt/spack/linux-ubuntu20.04-zen2/clang-13.0.1/cmake-3.22.2-r2xmbw62dyg2cidarrfycl7bjijwro55/bin/cmake

# The command to remove a file.
RM = /home/sajid/packages/spack/opt/spack/linux-ubuntu20.04-zen2/clang-13.0.1/cmake-3.22.2-r2xmbw62dyg2cidarrfycl7bjijwro55/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/sajid/packages/minigia

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/sajid/packages/minigia/build

# Include any dependencies generated for this target.
include src/utils/CMakeFiles/synergia_serialization.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/utils/CMakeFiles/synergia_serialization.dir/compiler_depend.make

# Include the progress variables for this target.
include src/utils/CMakeFiles/synergia_serialization.dir/progress.make

# Include the compile flags for this target's objects.
include src/utils/CMakeFiles/synergia_serialization.dir/flags.make

src/utils/CMakeFiles/synergia_serialization.dir/cereal_files.cc.o: src/utils/CMakeFiles/synergia_serialization.dir/flags.make
src/utils/CMakeFiles/synergia_serialization.dir/cereal_files.cc.o: ../src/utils/cereal_files.cc
src/utils/CMakeFiles/synergia_serialization.dir/cereal_files.cc.o: src/utils/CMakeFiles/synergia_serialization.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sajid/packages/minigia/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/utils/CMakeFiles/synergia_serialization.dir/cereal_files.cc.o"
	cd /home/sajid/packages/minigia/build/src/utils && /home/sajid/packages/spack/opt/spack/linux-ubuntu20.04-zen2/gcc-11.2.0/llvm-13.0.1-mpuqj67i6qekqo7eo74no2ahuvx4idli/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/utils/CMakeFiles/synergia_serialization.dir/cereal_files.cc.o -MF CMakeFiles/synergia_serialization.dir/cereal_files.cc.o.d -o CMakeFiles/synergia_serialization.dir/cereal_files.cc.o -c /home/sajid/packages/minigia/src/utils/cereal_files.cc

src/utils/CMakeFiles/synergia_serialization.dir/cereal_files.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/synergia_serialization.dir/cereal_files.cc.i"
	cd /home/sajid/packages/minigia/build/src/utils && /home/sajid/packages/spack/opt/spack/linux-ubuntu20.04-zen2/gcc-11.2.0/llvm-13.0.1-mpuqj67i6qekqo7eo74no2ahuvx4idli/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sajid/packages/minigia/src/utils/cereal_files.cc > CMakeFiles/synergia_serialization.dir/cereal_files.cc.i

src/utils/CMakeFiles/synergia_serialization.dir/cereal_files.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/synergia_serialization.dir/cereal_files.cc.s"
	cd /home/sajid/packages/minigia/build/src/utils && /home/sajid/packages/spack/opt/spack/linux-ubuntu20.04-zen2/gcc-11.2.0/llvm-13.0.1-mpuqj67i6qekqo7eo74no2ahuvx4idli/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sajid/packages/minigia/src/utils/cereal_files.cc -o CMakeFiles/synergia_serialization.dir/cereal_files.cc.s

# Object files for target synergia_serialization
synergia_serialization_OBJECTS = \
"CMakeFiles/synergia_serialization.dir/cereal_files.cc.o"

# External object files for target synergia_serialization
synergia_serialization_EXTERNAL_OBJECTS =

src/utils/libsynergia_serialization.a: src/utils/CMakeFiles/synergia_serialization.dir/cereal_files.cc.o
src/utils/libsynergia_serialization.a: src/utils/CMakeFiles/synergia_serialization.dir/build.make
src/utils/libsynergia_serialization.a: src/utils/CMakeFiles/synergia_serialization.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/sajid/packages/minigia/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libsynergia_serialization.a"
	cd /home/sajid/packages/minigia/build/src/utils && $(CMAKE_COMMAND) -P CMakeFiles/synergia_serialization.dir/cmake_clean_target.cmake
	cd /home/sajid/packages/minigia/build/src/utils && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/synergia_serialization.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/utils/CMakeFiles/synergia_serialization.dir/build: src/utils/libsynergia_serialization.a
.PHONY : src/utils/CMakeFiles/synergia_serialization.dir/build

src/utils/CMakeFiles/synergia_serialization.dir/clean:
	cd /home/sajid/packages/minigia/build/src/utils && $(CMAKE_COMMAND) -P CMakeFiles/synergia_serialization.dir/cmake_clean.cmake
.PHONY : src/utils/CMakeFiles/synergia_serialization.dir/clean

src/utils/CMakeFiles/synergia_serialization.dir/depend:
	cd /home/sajid/packages/minigia/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/sajid/packages/minigia /home/sajid/packages/minigia/src/utils /home/sajid/packages/minigia/build /home/sajid/packages/minigia/build/src/utils /home/sajid/packages/minigia/build/src/utils/CMakeFiles/synergia_serialization.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/utils/CMakeFiles/synergia_serialization.dir/depend
