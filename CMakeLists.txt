# CMakeList.txt : CMake project for FlashSimulations, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)
project(Flash_Simulation CXX)

add_subdirectory("Utilities")
add_subdirectory("Ideal Solutions")
add_subdirectory("Non Ideal Solutions")
add_subdirectory("Differential Vaporization")
add_subdirectory("Fugacity Pure Substances")
add_subdirectory("Mixture Equilibrium EOS")

