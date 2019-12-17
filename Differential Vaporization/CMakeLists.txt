# CMakeList.txt : CMake project for FlashSimulations, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

add_executable(DiffVap DiffVap.cpp)
add_executable(DiffVapEOS DiffVapEOS.cpp)


target_link_libraries(DiffVap Utilities)
target_link_libraries(DiffVapEOS Utilities)

