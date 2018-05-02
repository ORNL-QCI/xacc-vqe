#***********************************************************************************
# Copyright (c) 2016, UT-Battelle
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#   * Neither the name of the xacc nor the
#     names of its contributors may be used to endorse or promote products
#     derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contributors:
#   Initial API and implementation - Alex McCaskey
#
#**********************************************************************************/

set (PACKAGE_NAME "XACC VQE Runtime")
set (PACKAGE_DESCIPTION "XACC Variational Quantum Eigensolver IR Support")

set (IR_LIBRARY_NAME xacc-vqe-ir)

file (GLOB_RECURSE HEADERS *.hpp)
file (GLOB SRC *.cpp algorithms/uccsd/*.cpp)

# Set up dependencies to resources to track changes
usFunctionGetResourceSource(TARGET ${IR_LIBRARY_NAME} OUT SRC)
# Generate bundle initialization code
usFunctionGenerateBundleInit(TARGET ${IR_LIBRARY_NAME} OUT SRC)

add_library(${IR_LIBRARY_NAME} SHARED ${SRC})

set(_bundle_name xacc_vqe_ir)

set_target_properties(${IR_LIBRARY_NAME} PROPERTIES
  # This is required for every bundle
  COMPILE_DEFINITIONS US_BUNDLE_NAME=${_bundle_name}
  # This is for convenience, used by other CMake functions
  US_BUNDLE_NAME ${_bundle_name}
  )

if(BUILD_SHARED_LIBS)
  set_target_properties(${IR_LIBRARY_NAME} PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1
    )
endif()

# Embed meta-data from a manifest.json file
usFunctionEmbedResources(TARGET ${IR_LIBRARY_NAME}
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_SOURCE_DIR}
  FILES
    manifest.json
  )

target_link_libraries(${IR_LIBRARY_NAME} ${XACC_LIBRARIES})

if(APPLE)
	set_target_properties(${IR_LIBRARY_NAME} PROPERTIES INSTALL_RPATH "@loader_path/../lib;@loader_path")
	set_target_properties(${IR_LIBRARY_NAME} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
else()
	set_target_properties(${IR_LIBRARY_NAME} PROPERTIES INSTALL_RPATH "$ORIGIN/../lib:$ORIGIN")
	set_target_properties(${IR_LIBRARY_NAME} PROPERTIES LINK_FLAGS "-shared")
endif()

install(TARGETS ${IR_LIBRARY_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/plugins)
install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/vqe)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/algorithms/uccsd)
include_directories(${CMAKE_SOURCE_DIR}/vqe/transformation/jw)

# Gather tests
file (GLOB test_files tests/*.cpp)
add_tests("${test_files}" "${CMAKE_CURRENT_SOURCE_DIR}" "${XACC_LIBRARIES};xacc-vqe-ir;pthread;dl")

