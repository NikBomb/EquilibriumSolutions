# CMakeList.txt : CMake project for FlashSimulations, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)


add_library(Utilities 
   src/bisection.cpp
   include/gases.hpp
   include/idealSolutions.hpp
   include/nonIdealSolutions.hpp)


target_include_directories(Utilities PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:include/> 
)