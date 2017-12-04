#
# Module that checks whether ARPACK++ is available and usable.
#
# Variables used by this module which you may want to set:
# ARPACKPP_ROOT           Path list to search for ARPACK++.
#
# Sets the following variables:
# ARPACKPP_FOUND          True if ARPACK++ available.
# ARPACKPP_INCLUDE_DIRS   Path to the ARPACK++ include directories.
# ARPACKPP_LIBRARIES      Link against these libraries to use ARPACK++.
#

# find ARPACK which is required by ARPACK++
find_package(ARPACK)
if(NOT(ARPACK_FOUND))
  message(WARNING "ARPACK not found, skipping search for ARPACK++.")
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determing location of ARPACK++ failed:\n"
    "ARPACK is required but could not be found\n\n")
  return()
endif(NOT(ARPACK_FOUND))

# look for header files, only at positions given by the user
find_path(ARPACKPP_INCLUDE_DIR
  NAMES "arssym.h"
  PATHS ${ARPACKPP_PREFIX} ${ARPACKPP_ROOT}
  PATH_SUFFIXES "include" "include/arpack++"
  NO_DEFAULT_PATH
)

# look for header files, including default paths
find_path(ARPACKPP_INCLUDE_DIR
  NAMES "arssym.h"
  PATH_SUFFIXES "include" "include/arpack++"
)

# check header usability
include(CMakePushCheckState)
cmake_push_check_state()

# we need if clauses here because variable is set variable-NOTFOUND if the
# searches above were not successful; without them CMake print errors like:
# "CMake Error: The following variables are used in this project, but they
# are set to NOTFOUND. Please set them or make sure they are set and tested
# correctly in the CMake files."
if(ARPACKPP_INCLUDE_DIR)
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${ARPACKPP_INCLUDE_DIR})
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${ARPACK_LIBRARIES}) #${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${UMFPACK_LIBRARIES} ${SUPERLU_LIBRARIES}
endif(ARPACKPP_INCLUDE_DIR)

# end of header usability check
cmake_pop_check_state()

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "ARPACKPP"
  DEFAULT_MSG
  ARPACKPP_INCLUDE_DIR
)

# hide the introduced cmake cached variables in cmake GUIs
mark_as_advanced(ARPACKPP_INCLUDE_DIR)

# if headers are found, store results
if(ARPACKPP_FOUND)
  set(ARPACKPP_INCLUDE_DIRS ${ARPACKPP_INCLUDE_DIR})
  set(ARPACKPP_LIBRARIES ${ARPACK_LIBRARIES}) #${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${UMFPACK_LIBRARIES} ${SUPERLU_LIBRARIES}
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determing location of ARPACK++ succeded:\n"
    "Include directory: ${ARPACKPP_INCLUDE_DIRS}\n"
    "Libraries to link against: ${ARPACKPP_LIBRARIES}\n\n")
  set(ARPACKPP_DUNE_COMPILE_FLAGS "-I${ARPACKPP_INCLUDE_DIRS}"
    CACHE STRING "Compile flags used by DUNE when compiling ARPACK++ programs")
  set(ARPACKPP_DUNE_LIBRARIES ${ARPACKPP_LIBRARIES}
    CACHE STRING "Libraries used by DUNE when linking ARPACK++ programs")
else(ARPACKPP_FOUND)
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determing location of ARPACK++ failed:\n"
    "Include directory: ${ARPACKPP_INCLUDE_DIRS}\n"
    "Libraries to link against: ${ARPACKPP_LIBRARIES}\n\n")
endif(ARPACKPP_FOUND)

# set HAVE_ARPACKPP for config.h
set(HAVE_ARPACKPP ${ARPACKPP_FOUND})

