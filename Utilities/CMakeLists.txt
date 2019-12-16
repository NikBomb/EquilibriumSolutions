# CMakeList.txt : CMake project for FlashSimulations, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)


add_library(Utilities
   include/gases.hpp
   include/idealSolutions.hpp
   include/nonIdealSolutions.hpp
   src/bisection.cpp
   src/poly34.cpp
   src/eos.cpp)


target_include_directories(Utilities PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:include/> 
)