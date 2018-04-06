# - Config file for XACC
# To point to your install of XACC, pass the 
# XACC_ROOT flag to your cmake configure.
#
# It defines the following variables
#  XACC_INCLUDE_DIRS - include directories for XACC
#  XACC_LIBRARIES    - libraries to link against
#  XACC_LIBRARY_DIR  - the XACC library directory 
if (NOT XACCVQE_ROOT)
   get_filename_component(XACCVQE_ROOT "${CMAKE_CURRENT_LIST_FILE}" PATH)
endif()
set (XACCVQE_INCLUDE_DIRS "${XACCVQE_ROOT}/include/vqe;${XACCVQE_ROOT}/include/cppoptlib")
set (XACCVQE_LIBRARY_DIR "${XACCVQE_ROOT}/lib;${XACCVQE_ROOT}/plugins/ir")
set(CppUsLib CppMicroServicesd)
link_directories("${XACCVQE_ROOT}/lib;${XACCVQE_ROOT}/plugins/ir;${XACCVQE_ROOT}/plugins/misc")
set (XACCVQE_LIBRARIES "xacc-vqe-ir;xacc-vqe-tasks")
