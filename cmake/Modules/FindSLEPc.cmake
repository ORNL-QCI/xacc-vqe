# - Try to find SLEPC
# Once done this will define
#
#  SLEPC_FOUND            - system has SLEPc
#  SLEPC_INCLUDE_DIRS     - include directories for SLEPc
#  SLEPC_LIBRARY_DIRS     - library directories for SLEPc
#  SLEPC_LIBARIES         - libraries for SLEPc
#  SLEPC_STATIC_LIBARIES  - ibraries for SLEPc (static linking, undefined if not required)
#  SLEPC_VERSION          - version of SLEPc
#  SLEPC_VERSION_MAJOR    - First number in SLEPC_VERSION
#  SLEPC_VERSION_MINOR    - Second number in SLEPC_VERSION
#  SLEPC_VERSION_SUBMINOR - Third number in SLEPC_VERSION

# Load CMake pkg-config module
find_package(PkgConfig REQUIRED)

# Find SLEPc pkg-config file
set(ENV{PKG_CONFIG_PATH} "$ENV{SLEPC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig:$ENV{SLEPC_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
set(ENV{PKG_CONFIG_PATH} "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig:$ENV{PETSC_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
set(ENV{PKG_CONFIG_PATH} "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}:$ENV{PETSC_DIR}:$ENV{PKG_CONFIG_PATH}")
pkg_search_module(SLEPC crayslepc_real SLEPc)

# Extract major, minor, etc from version string
if (SLEPC_VERSION)
  string(REPLACE "." ";" VERSION_LIST ${SLEPC_VERSION})
  list(GET VERSION_LIST 0 SLEPC_VERSION_MAJOR)
  list(GET VERSION_LIST 1 SLEPC_VERSION_MINOR)
  list(GET VERSION_LIST 2 SLEPC_VERSION_SUBMINOR)
endif()

# Configure SLEPc IMPORT (this involves creating an 'imported' target
# and attaching 'properties')
if (SLEPC_FOUND AND NOT TARGET SLEPC::slepc)
  add_library(SLEPC::slepc INTERFACE IMPORTED)

  # Add include paths
  set_property(TARGET SLEPC::slepc PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${SLEPC_INCLUDE_DIRS})

  # Add libraries
  unset(_libs)
  foreach (lib ${SLEPC_LIBRARIES})
    find_library(LIB_${lib} NAMES ${lib} PATHS ${SLEPC_LIBRARY_DIRS} NO_DEFAULT_PATH)
    list(APPEND _libs ${LIB_${lib}})
  endforeach()
  set_property(TARGET SLEPC::slepc PROPERTY INTERFACE_LINK_LIBRARIES "${_libs}")
endif()

if (SLEPC_FOUND AND NOT TARGET SLEPC::slepc_static)
  add_library(SLEPC::slepc_static INTERFACE IMPORTED)

  # Add libraries (static)
  unset(_libs)
  foreach (lib ${SLEPC_STATIC_LIBRARIES})
    find_library(LIB_${lib} ${lib} HINTS ${SLEPC_STATIC_LIBRARY_DIRS})
    list(APPEND _libs ${LIB_${lib}})
  endforeach()
  set_property(TARGET SLEPC::slepc_static PROPERTY INTERFACE_LINK_LIBRARIES "${_libs}")

endif()

# Standard package handling
include(FindPackageHandleStandardArgs)
if (SLEPC_FOUND)
  find_package_handle_standard_args(SLEPc
    REQUIRED_VARS SLEPC_FOUND 
    VERSION_VAR SLEPC_VERSION
    FAIL_MESSAGE "SLEPc could not be configured.")
else()
  find_package_handle_standard_args(SLEPc
    REQUIRED_VARS SLEPC_FOUND
    FAIL_MESSAGE "SLEPc could not be found. Be sure to set SLEPC_DIR.")
endif()