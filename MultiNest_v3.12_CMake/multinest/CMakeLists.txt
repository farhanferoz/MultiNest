cmake_minimum_required(VERSION 2.8)
project(multinest Fortran C CXX)

# By default CMake will install MultiNest to /usr/local on Linux/Mac
# unless you override the install location manually. You can do this
# either by either running:
#  cmake -DCMAKE_INSTALL_PREFIX=/path/to/dir ..
# or by uncommenting and modifying this variable:
#SET(CMAKE_INSTALL_PREFIX /usr/local/multinest)

#=============================================================================
# All done editing
SET(MultiNest_VERSION "3.11")

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_subdirectory(src)
