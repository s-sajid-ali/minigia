# Install script for directory: /home/sajid/packages/minigia/src/utils

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/home/sajid/packages/spack/opt/spack/linux-ubuntu20.04-zen2/gcc-11.2.0/llvm-13.0.1-mpuqj67i6qekqo7eo74no2ahuvx4idli/bin/llvm-objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsynergia_parallel_utils.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsynergia_parallel_utils.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsynergia_parallel_utils.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/sajid/packages/minigia/build/src/utils/libsynergia_parallel_utils.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsynergia_parallel_utils.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsynergia_parallel_utils.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsynergia_parallel_utils.so"
         OLD_RPATH "/home/sajid/packages/spack/var/spack/environments/synergia/.spack-env/view/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/home/sajid/packages/spack/opt/spack/linux-ubuntu20.04-zen2/gcc-11.2.0/llvm-13.0.1-mpuqj67i6qekqo7eo74no2ahuvx4idli/bin/llvm-strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libsynergia_parallel_utils.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/sajid/packages/minigia/build/src/utils/libsynergia_hdf5_utils.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/sajid/packages/minigia/build/src/utils/libsynergia_serialization.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/sajid/packages/minigia/build/src/utils/libsynergia_command_line.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/liblsexpr.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/liblsexpr.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/liblsexpr.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/sajid/packages/minigia/build/src/utils/liblsexpr.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/liblsexpr.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/liblsexpr.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/liblsexpr.so"
         OLD_RPATH "/home/sajid/packages/spack/var/spack/environments/synergia/.spack-env/view/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/home/sajid/packages/spack/opt/spack/linux-ubuntu20.04-zen2/gcc-11.2.0/llvm-13.0.1-mpuqj67i6qekqo7eo74no2ahuvx4idli/bin/llvm-strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/liblsexpr.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/synergia/utils/comm_converter.h;/synergia/utils/command_line_arg.h;/synergia/utils/commxx.h;/synergia/utils/commxx_divider.h;/synergia/utils/fast_int_floor.h;/synergia/utils/floating_point.h;/synergia/utils/gsvector.h;/synergia/utils/hdf5_misc.h;/synergia/utils/hdf5_file.h;/synergia/utils/hdf5_serial_writer.h;/synergia/utils/hdf5_writer.h;/synergia/utils/complex_error_function.h;/synergia/utils/multi_array_typedefs.h;/synergia/utils/parallel_utils.h;/synergia/utils/simple_timer.h;/synergia/utils/cereal.h;/synergia/utils/cereal_files.h;/synergia/utils/digits.h;/synergia/utils/logger.h;/synergia/utils/lsexpr.h;/synergia/utils/synergia_omp.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/synergia/utils" TYPE FILE FILES
    "/home/sajid/packages/minigia/src/utils/comm_converter.h"
    "/home/sajid/packages/minigia/src/utils/command_line_arg.h"
    "/home/sajid/packages/minigia/src/utils/commxx.h"
    "/home/sajid/packages/minigia/src/utils/commxx_divider.h"
    "/home/sajid/packages/minigia/src/utils/fast_int_floor.h"
    "/home/sajid/packages/minigia/src/utils/floating_point.h"
    "/home/sajid/packages/minigia/src/utils/gsvector.h"
    "/home/sajid/packages/minigia/src/utils/hdf5_misc.h"
    "/home/sajid/packages/minigia/src/utils/hdf5_file.h"
    "/home/sajid/packages/minigia/src/utils/hdf5_serial_writer.h"
    "/home/sajid/packages/minigia/src/utils/hdf5_writer.h"
    "/home/sajid/packages/minigia/src/utils/complex_error_function.h"
    "/home/sajid/packages/minigia/src/utils/multi_array_typedefs.h"
    "/home/sajid/packages/minigia/src/utils/parallel_utils.h"
    "/home/sajid/packages/minigia/src/utils/simple_timer.h"
    "/home/sajid/packages/minigia/src/utils/cereal.h"
    "/home/sajid/packages/minigia/src/utils/cereal_files.h"
    "/home/sajid/packages/minigia/src/utils/digits.h"
    "/home/sajid/packages/minigia/src/utils/logger.h"
    "/home/sajid/packages/minigia/src/utils/lsexpr.h"
    "/home/sajid/packages/minigia/src/utils/synergia_omp.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/sajid/packages/minigia/build/src/utils/tests/cmake_install.cmake")
endif()

