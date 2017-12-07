
find_package(PETSc) 
find_package(SLEPc)

if (PETSC_FOUND AND SLEPC_FOUND) 

set (PACKAGE_NAME "XACC Slepc Ground State Energy Calculator")
set (PACKAGE_DESCIPTION "")

set (LIBRARY_NAME xacc-vqe-slepc)

set(CMAKE_INSTALL_RPATH "${XACC_ROOT}/lib;${XACC_ROOT}/lib/plugins/misc;/usr/local/lib64")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/task/tasks)

include_directories(/projects/petsc/arch-linux2-c-debug/include)
include_directories(/projects/petsc/include)
include_directories(/projects/slepc-3.8.2/arch-linux2-c-debug/include)
include_directories(/projects/slepc-3.8.2/include)

file (GLOB_RECURSE HEADERS *.hpp)

file (GLOB SRC *.cpp)

# Set up dependencies to resources to track changes
usFunctionGetResourceSource(TARGET ${LIBRARY_NAME} OUT SRC)
# Generate bundle initialization code
usFunctionGenerateBundleInit(TARGET ${LIBRARY_NAME} OUT SRC)

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

target_link_libraries(${LIBRARY_NAME} ${Boost_LIBRARIES} ${XACC_LIBRARIES} ${PYTHON_LIBRARY} xacc-vqe-ir ${MPI_CXX_LIBRARIES} /projects/petsc/arch-linux2-c-debug/lib/libpetsc.so /projects/slepc-3.8.2/arch-linux2-c-debug/lib/libslepc.so)

install(TARGETS ${LIBRARY_NAME} DESTINATION ${XACC_LIBRARY_DIR}/plugins/misc)

# Gather tests
#file (GLOB test_files tests/*.cpp)
#add_tests("${test_files}" "${CMAKE_CURRENT_SOURCE_DIR}" "${Boost_LIBRARIES};pthread;dl;xacc-vqe-ir;xacc-vqe-tasks;xacc-vqe-irtransformations")

endif()