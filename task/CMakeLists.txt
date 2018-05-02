set (PACKAGE_NAME "XACC VQE Tasks")
set (PACKAGE_DESCIPTION "XACC Variational Quantum Eigensolver Tasks")

set (LIBRARY_NAME xacc-vqe-tasks)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/tasks)

file (GLOB_RECURSE HEADERS *.hpp)

file (GLOB SRC *.cpp tasks/*.cpp)

# Set up dependencies to resources to track changes
usFunctionGetResourceSource(TARGET ${LIBRARY_NAME} OUT SRC)
# Generate bundle initialization code
usFunctionGenerateBundleInit(TARGET ${LIBRARY_NAME} OUT SRC)

add_library(${LIBRARY_NAME} SHARED ${SRC})

set(_bundle_name xacc_vqe_tasks)

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

target_link_libraries(${LIBRARY_NAME} ${XACC_LIBRARIES} xacc-vqe-ir)

if(APPLE)
   set_target_properties(${LIBRARY_NAME} PROPERTIES INSTALL_RPATH "@loader_path/../lib;@loader_path")
   set_target_properties(${LIBRARY_NAME} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
else()
   set_target_properties(${LIBRARY_NAME} PROPERTIES INSTALL_RPATH "$ORIGIN/../lib:$ORIGIN")
   set_target_properties(${LIBRARY_NAME} PROPERTIES LINK_FLAGS "-shared")
endif()

install(TARGETS ${LIBRARY_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/plugins)
install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/vqe)

# Gather tests
file (GLOB test_files tests/*.cpp)
add_tests("${test_files}" "${CMAKE_CURRENT_SOURCE_DIR}" "${XACC_LIBRARIES};pthread;dl;xacc-vqe-ir;xacc-vqe-tasks;xacc-vqe-irtransformations")

add_subdirectory(tasks)
