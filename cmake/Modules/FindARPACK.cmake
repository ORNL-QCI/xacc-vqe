#
# Module that checks whether ARPACK is available and usable.
#
# Variables used by this module which you may want to set:
# ARPACK_ROOT           Path list to search for ARPACK.
#
# Sets the following variables:
# ARPACK_FOUND          True if ARPACK available.
# ARPACK_LIBRARIES      Link against these libraries to use ARPACK.
#

# check for Fortran support which is required by ARPACK
#if(NOT(Fortran_Works))
#  message(WARNING "Fortran doesn't seem to be supported, skipping search for ARPACK.")
  # log errornous result
#  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
#    "Determing location of ARPACK failed:\n"
#    "Fortran support is required but could not be detected\n\n")
#  return()
#endif(NOT(Fortran_Works))

# look for library, only at positions given by the user
find_library(ARPACK_LIBRARY
  NAMES "arpack"
  PATHS ${ARPACK_PREFIX} ${ARPACK_ROOT}
  PATH_SUFFIXES "lib" "lib32" "lib64"
  NO_DEFAULT_PATH
)

# look for library files, including default paths
find_library(ARPACK_LIBRARY
  NAMES "arpack"
  PATH_SUFFIXES "lib" "lib32" "lib64"
)

# check header usability
include(CMakePushCheckState)
cmake_push_check_state()

# we need if clauses here because variable is set variable-NOTFOUND if the
# searches above were not successful; without them CMake print errors like:
# "CMake Error: The following variables are used in this project, but they
# are set to NOTFOUND. Please set them or make sure they are set and tested
# correctly in the CMake files."
if(ARPACK_LIBRARY)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${ARPACK_LIBRARY})
endif(ARPACK_LIBRARY)

# end of header usability check
cmake_pop_check_state()

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "ARPACK"
  DEFAULT_MSG
  ARPACK_LIBRARY
)

# hide the introduced cmake cached variables in cmake GUIs
mark_as_advanced(ARPACK_LIBRARY)

# if headers are found, store results
if(ARPACK_FOUND)
  set(ARPACK_LIBRARIES ${ARPACK_LIBRARY})
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determing location of ARPACK succeded:\n"
    "Libraries to link against: ${ARPACK_LIBRARIES}\n\n")
  set(ARPACK_DUNE_LIBRARIES ${ARPACK_LIBRARIES}
    CACHE STRING "Libraries used by DUNE when linking ARPACK programs")
else(ARPACK_FOUND)
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determing location of ARPACK failed:\n"
    "Libraries to link against: ${ARPACK_LIBRARIES}\n\n")
endif(ARPACK_FOUND)
