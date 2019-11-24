cmake_minimum_required(VERSION 2.8) 

#=============================================================================
# Copyright 2012 Brian Kloppenborg
#
# This code is licensed under the MIT License.  
#
# The MIT License
#
# License for the specific language governing rights and limitations under
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
#=============================================================================

# Find required packages
FIND_PACKAGE(LAPACK REQUIRED)

# If we are compiling with gfortran, add -ffree-line-length-none
foreach(lang C CXX Fortran)
    if(CMAKE_${lang}_COMPILER_ID STREQUAL GNU)
        MESSAGE(STATUS "Detected gfortran, adding -ffree-line-length-none compiler flag.")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
        break()
    endif()
endforeach()

# Emit a warning for gfortran < 4.3.6 and cwrapper.f90
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
 if(${CMAKE_C_COMPILER_VERSION} VERSION_LESS 4.3.7)
    message(WARNING "Insufficient gfortran version for the cwrapper module. PyMultiNest will not work automagically. Please refer to the PyMultiNest manual if you want to use it. https://johannesbuchner.github.io/PyMultiNest/install.html")
 endif()
endif()

# Enable optimization.
#LIST(APPEND CMAKE_Fortran_FLAGS "-O3")

# ============================
# Collect the source:
# ============================
FILE(GLOB MultiNest_SOURCE *.f90 *.F90)
# In previous versions of MultiNest (<3.8) the cwrapper.f90 file presented
# some issues with older compilers, thus it was disabled by default. This
# was fixed in r34, thus it is enabled by default now. Enable the following
# two lines if it presents you with trouble.
#FILE(GLOB MultiNest_REMOVE cwrapper.f90)
#list(REMOVE_ITEM MultiNest_SOURCE ${MultiNest_REMOVE})

# ============================
# Build the non-MPI libraries:
# ============================
# libmultinest.a
ADD_LIBRARY(multinest_static STATIC ${MultiNest_SOURCE})
SET_TARGET_PROPERTIES(multinest_static PROPERTIES LINKER_LANGUAGE Fortran 
    VERSION ${MultiNest_VERSION} OUTPUT_NAME multinest
)
TARGET_LINK_LIBRARIES(multinest_static ${LAPACK_LIBRARIES})

# libmultinest.so
ADD_LIBRARY(multinest_shared SHARED ${MultiNest_SOURCE})
SET_TARGET_PROPERTIES(multinest_shared PROPERTIES LINKER_LANGUAGE Fortran 
    VERSION ${MultiNest_VERSION} OUTPUT_NAME multinest
)
TARGET_LINK_LIBRARIES(multinest_shared ${LAPACK_LIBRARIES})  

# Set the default library against which to link:
SET(MultiNest_LIBRARIES multinest_shared)

# ============================
# Build MPI-enabled libraries, if MPI is installed
# ============================
FIND_PACKAGE(MPI)
if(NOT MPI_Fortran_FOUND)
    MESSAGE(STATUS "MPI not found, only non-MPI MultiNest libraries will be built.")
else()
    ADD_LIBRARY(multinest_mpi_static STATIC ${MultiNest_SOURCE})
    SET_TARGET_PROPERTIES(multinest_mpi_static PROPERTIES LINKER_LANGUAGE Fortran VERSION ${MultiNest_VERSION} OUTPUT_NAME multinest_mpi)
    TARGET_LINK_LIBRARIES(multinest_mpi_static ${LAPACK_LIBRARIES} ${MPI_Fortran_LIBRARIES})
    
    ADD_LIBRARY(multinest_mpi_shared SHARED ${MultiNest_SOURCE})
    SET_TARGET_PROPERTIES(multinest_mpi_shared PROPERTIES LINKER_LANGUAGE Fortran VERSION ${MultiNest_VERSION} OUTPUT_NAME multinest_mpi)    
    TARGET_LINK_LIBRARIES(multinest_mpi_shared ${LAPACK_LIBRARIES} ${MPI_Fortran_LIBRARIES})
    
    # Include directories, set linking requirements
    INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})

    # Set compilation flags, ensure -DMPI is, at a minimum, set.
    if(${MPI_Fortran_COMPILE_FLAGS})
        SET_TARGET_PROPERTIES(multinest_mpi_static multinest_mpi_shared PROPERTIES COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS}) 
    else()
        SET_TARGET_PROPERTIES(multinest_mpi_static multinest_mpi_shared PROPERTIES COMPILE_FLAGS "-DMPI")
    endif()
    
    if(${MPI_Fortran_LINK_FLAGS})
        SET_TARGET_PROPERTIES(multinest_mpi_static multinest_mpi_shared PROPERTIES LINK_FLAGS ${MPI_Fortran_LINK_FLAGS}) 
    endif()
    
    # Set the default linking option to use the shared library.
    SET(MultiNest_LIBRARIES multinest_mpi_shared)
endif ()

# ============================
# After compilation, copy the Fortran modules into the 'modules' directory.
# ============================
ADD_CUSTOM_COMMAND(
    TARGET multinest_static
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/src/kmeans_clstr.mod ${CMAKE_SOURCE_DIR}/modules/
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/src/nested.mod ${CMAKE_SOURCE_DIR}/modules/
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/src/posterior.mod ${CMAKE_SOURCE_DIR}/modules/
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/src/priors.mod ${CMAKE_SOURCE_DIR}/modules/
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/src/randomns.mod ${CMAKE_SOURCE_DIR}/modules/
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/src/utils1.mod ${CMAKE_SOURCE_DIR}/modules/
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/src/xmeans_clstr.mod ${CMAKE_SOURCE_DIR}/modules/
# cnested is not built by default anymore so it need not be copied
#    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/src/cnested.mod ${CMAKE_SOURCE_DIR}/modules/
)

# ============================
# Installation commands
# ============================

# We need to modify the RPATH between the compilation and installation phases.
# When we first build the programs, we use the full RPATH within the build tree
# instead of the default system RPATH. This ensures the binaries are executable
# right away.
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# When we install, switch the RPATH to the default system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

# Determine which libraries we will be installing. 
SET(MultiNest_INSTALL_LIBS multinest_static multinest_shared)
if(MPI_Fortran_FOUND)
    LIST(APPEND MultiNest_INSTALL_LIBS multinest_mpi_static multinest_mpi_shared)
endif(MPI_Fortran_FOUND)

# Specify where compiled files should go
INSTALL(TARGETS ${MultiNest_INSTALL_LIBS} DESTINATION lib)
# Don't forget the fortran modules!
INSTALL(FILES 
    ${CMAKE_BINARY_DIR}/src/kmeans_clstr.mod
    ${CMAKE_BINARY_DIR}/src/nested.mod
    ${CMAKE_BINARY_DIR}/src/posterior.mod
    ${CMAKE_BINARY_DIR}/src/priors.mod 
    ${CMAKE_BINARY_DIR}/src/randomns.mod
    ${CMAKE_BINARY_DIR}/src/utils1.mod
    ${CMAKE_BINARY_DIR}/src/xmeans_clstr.mod
# cnested is not built by default anymore so it need not be installed.
#    ${CMAKE_BINARY_DIR}/src/cnested.mod
    DESTINATION modules
)
# Install the multinest.h file
# Install the FindMultiNest.cmake script (this probably isn't the right place)
INSTALL(FILES ${CMAKE_SOURCE_DIR}/include/multinest.h DESTINATION include)
INSTALL(FILES ${CMAKE_SOURCE_DIR}/CMakeModules/FindMultiNest.cmake DESTINATION
    CMakeModules
)
# and ensure the chains directory is created too
INSTALL(DIRECTORY
    ${CMAKE_SOURCE_DIR}/bin/chains
    DESTINATION bin
)

# ============================
# Now with the libraries specified, build all of the executables:
# ============================
ADD_SUBDIRECTORY(example_ackley)
ADD_SUBDIRECTORY(example_eggbox_C)
ADD_SUBDIRECTORY(example_eggbox_C++)
ADD_SUBDIRECTORY(example_gaussian)
ADD_SUBDIRECTORY(example_gauss_shell)
ADD_SUBDIRECTORY(example_himmelblau)
ADD_SUBDIRECTORY(example_obj_detect)
ADD_SUBDIRECTORY(example_rosenbrock)
