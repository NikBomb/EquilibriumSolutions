# CMakeList.txt : CMake project for FlashSimulations, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

add_executable(NonIdeal NonIdeal.cpp)

target_link_libraries(NonIdeal Utilities)

