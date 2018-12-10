set (PACKAGE_NAME "XACC Slepc Ground State Energy Calculator")
set (PACKAGE_DESCIPTION "")

set (LIBRARY_NAME xacc-vqe-slepc)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/task/tasks)

message(STATUS "${PETSC_INCLUDE_DIRS}")
include_directories(${MPI_CXX_INCLUDE_PATH})
include_directories(${PETSC_INCLUDE_DIRS})
include_directories(${SLEPC_INCLUDE_DIRS})

file (GLOB_RECURSE HEADERS *.hpp)

file (GLOB SRC *.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPETSC_USE_COMPLEX")
# Set up dependencies to resources to track changes
usFunctionGetResourceSource(TARGET ${LIBRARY_NAME} OUT SRC)
# Generate bundle initialization code
usFunctionGenerateBundleInit(TARGET ${LIBRARY_NAME} OUT SRC)

link_directories(${PETSC_LIBRARY_DIRS})
link_directories(${SLEPC_LIBRARY_DIRS})
add_library(${LIBRARY_NAME} SHARED ${SRC})

set(_bundle_name xacc_vqe_slepc)

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

target_link_libraries(${LIBRARY_NAME} ${XACC_LIBRARIES} ${PYTHON_LIBRARY} xacc-vqe-ir ${MPI_CXX_LIBRARIES} ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES})
if(APPLE)
   set_target_properties(${LIBRARY_NAME} PROPERTIES INSTALL_RPATH "@loader_path/../lib;@loader_path")
   set_target_properties(${LIBRARY_NAME} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
else()
   set_target_properties(${LIBRARY_NAME} PROPERTIES INSTALL_RPATH "$ORIGIN/../lib:$ORIGIN")
   set_target_properties(${LIBRARY_NAME} PROPERTIES LINK_FLAGS "-shared")
endif()

install(TARGETS ${LIBRARY_NAME} DESTINATION $ENV{HOME}/.xacc/plugins)

# Gather tests
#file (GLOB test_files tests/*.cpp)
#add_tests("${test_files}" "${CMAKE_CURRENT_SOURCE_DIR}" "${Boost_LIBRARIES};pthread;dl;xacc-vqe-ir;xacc-vqe-tasks;xacc-vqe-irtransformations")

