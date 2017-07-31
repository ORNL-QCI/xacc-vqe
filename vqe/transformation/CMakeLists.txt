set (PACKAGE_NAME "XACC VQE IR Transformations")
set (PACKAGE_DESCIPTION "XACC Variational Quantum Eigensolver IR Transformations")

set (LIBRARY_NAME xacc-vqe-irtransformations)

file (GLOB_RECURSE HEADERS *.hpp)
file (GLOB SRC *.cpp)

# Set up dependencies to resources to track changes
usFunctionGetResourceSource(TARGET ${LIBRARY_NAME} OUT SRC)
# Generate bundle initialization code
usFunctionGenerateBundleInit(TARGET ${LIBRARY_NAME} OUT SRC)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(${LIBRARY_NAME} SHARED ${SRC})

set(_bundle_name xacc_vqe_irtransformations)

set_target_properties(${LIBRARY_NAME} PROPERTIES
  # This is required for every bundle
  COMPILE_DEFINITIONS US_BUNDLE_NAME=${_bundle_name}
  # This is for convenience, used by other CMake functions
  US_BUNDLE_NAME ${_bundle_name}
  )

if(BUILD_SHARED_LIBS)
  set_target_properties(${LIBRARY_NAME} PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1
    )
endif()

# Embed meta-data from a manifest.json file
usFunctionEmbedResources(TARGET ${LIBRARY_NAME}
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_SOURCE_DIR}
  FILES
    manifest.json
  )

target_link_libraries(${LIBRARY_NAME} ${Boost_LIBRARIES} ${XACC_LIBRARIES} ${SYMENGINE_LIBRARIES})

install(TARGETS ${LIBRARY_NAME} DESTINATION ${XACC_LIBRARY_DIR}/plugins/transformations)

# Gather tests
file (GLOB test_files tests/*.cpp)
add_tests("${test_files}" "${CMAKE_CURRENT_SOURCE_DIR}" "${Boost_LIBRARIES};pthread;dl;xacc-vqe-ir;xacc-vqe-irtransformations")

