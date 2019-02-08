# Install script for directory: /home/qlu/s2-devel/build/synergia2/src/synergia/utils/vectorclass

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/qlu/s2-devel/install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
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
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/synergia/utils/vectorclass" TYPE FILE FILES
    "/home/qlu/s2-devel/build/synergia2/src/synergia/utils/vectorclass/instrset.h"
    "/home/qlu/s2-devel/build/synergia2/src/synergia/utils/vectorclass/instrset_detect.cpp"
    "/home/qlu/s2-devel/build/synergia2/src/synergia/utils/vectorclass/license.txt"
    "/home/qlu/s2-devel/build/synergia2/src/synergia/utils/vectorclass/vectorclass.h"
    "/home/qlu/s2-devel/build/synergia2/src/synergia/utils/vectorclass/vectorf128.h"
    "/home/qlu/s2-devel/build/synergia2/src/synergia/utils/vectorclass/vectorf256.h"
    "/home/qlu/s2-devel/build/synergia2/src/synergia/utils/vectorclass/vectorf256e.h"
    "/home/qlu/s2-devel/build/synergia2/src/synergia/utils/vectorclass/vectorf512.h"
    "/home/qlu/s2-devel/build/synergia2/src/synergia/utils/vectorclass/vectorf512e.h"
    "/home/qlu/s2-devel/build/synergia2/src/synergia/utils/vectorclass/vectori128.h"
    "/home/qlu/s2-devel/build/synergia2/src/synergia/utils/vectorclass/vectori256.h"
    "/home/qlu/s2-devel/build/synergia2/src/synergia/utils/vectorclass/vectori256e.h"
    "/home/qlu/s2-devel/build/synergia2/src/synergia/utils/vectorclass/vectori512.h"
    "/home/qlu/s2-devel/build/synergia2/src/synergia/utils/vectorclass/vectori512e.h"
    "/home/qlu/s2-devel/build/synergia2/src/synergia/utils/vectorclass/vectormath_common.h"
    "/home/qlu/s2-devel/build/synergia2/src/synergia/utils/vectorclass/vectormath_exp.h"
    "/home/qlu/s2-devel/build/synergia2/src/synergia/utils/vectorclass/vectormath_hyp.h"
    "/home/qlu/s2-devel/build/synergia2/src/synergia/utils/vectorclass/vectormath_lib.h"
    "/home/qlu/s2-devel/build/synergia2/src/synergia/utils/vectorclass/vectormath_trig.h"
    )
endif()

