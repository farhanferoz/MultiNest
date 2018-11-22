# - Find MultiNest

# Defines the following variables:
# MULTINEST_INCLUDE_DIRS    - Location of MULTINEST's include directory.
# MULTINEST_LIBRARIES       - Location of MULTINEST's libraries
# MULTINEST_FOUND           - True if MULTINEST has been located
# MULTINEST_MODULES         - Location of MULTINEST's Fortran modules
#
# The following components may be requested:
#   MPI - MULTINEST library with MPI enabled.
#
# You may provide a hint to where MULTINEST's root directory may be located
# by setting MULTINEST_ROOT before calling this script.
#
# Variables used by this module, they can change the default behaviour and
# need to be set before calling find_package:
#
#=============================================================================
#
#   MULTINEST_USE_MPI   Can be set to ON/true to find the MPI-enabled version of
#                       MULTINEST.  Defaults to false.
#
#=============================================================================
# Copyright 2012 Brian Kloppenborg
#
#  This code is licensed under the MIT License.  See the FindMULTINEST.cmake script
#  for the text of the license.
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

FIND_PACKAGE(LAPACK REQUIRED)
SET(_INCLUDES ${LAPACK_INCLUDES})
SET(_LIBS ${LAPACK_LIBRARIES})

IF(MULTINEST_INCLUDES)
  # Already in cache, be silent
  set (MULTINEST_FIND_QUIETLY TRUE)
ENDIF (MULTINEST_INCLUDES)

FIND_PATH(MULTINEST_ROOT_DIR
    NAMES include/multinest.h
    HINTS /usr/local/multinest ${MULTINEST_ROOT} 
    DOC "MULTINEST root directory.")
    
FIND_PATH(_MULTINEST_INCLUDE_DIRS
    NAMES multinest.h
    HINTS ${MULTINEST_ROOT_DIR}/include
    DOC "MULTINEST Include directory")
    
FIND_PATH(MULTINEST_MODULE_DIRS
    NAMES nested.mod
    HINTS /usr/local/modules ${MULTINEST_ROOT_DIR}/modules
    DOC "MULTINEST Fortran module directory")

# Now find the library:
SET(_MULTINEST_LIB_NAME "multinest")
IF(MULTINEST_USE_MPI)
    SET(_MULTINEST_LIB_NAME "multinest_mpi")
ENDIF(MULTINEST_USE_MPI)

FIND_LIBRARY(_MULTINEST_LIBRARY
    NAMES ${_MULTINEST_LIB_NAME}
    HINTS ${MULTINEST_ROOT_DIR}/lib)

SET(MULTINEST_INCLUDE_DIRS ${_MULTINEST_INCLUDE_DIRS} ${_INCLUDES})
SET(MULTINEST_LIBRARIES ${_MULTINEST_LIBRARY} ${_LIBS})

# handle the QUIETLY and REQUIRED arguments and set MULTINEST_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE (FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(MULTINEST DEFAULT_MSG MULTINEST_LIBRARIES MULTINEST_INCLUDE_DIRS MULTINEST_MODULE_DIRS)

MARK_AS_ADVANCED(MULTINEST_LIBRARIES MULTINEST_INCLUDE_DIRS MULTINEST_MODULE_DIRS)


