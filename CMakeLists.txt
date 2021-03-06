project(xacc-vqe CXX)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_DISABLE_IN_SOURCE_BUILDS ON)
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)

if (FROM_SETUP_PY AND NOT APPLE)
   message(STATUS "Running build from setup.py, linking to static libstdc++")
   set (CMAKE_SHARED_LINKER_FLAGS "-static-libstdc++" CACHE INTERNAL "" FORCE)
endif()
if(APPLE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -frtti -fexceptions -Wno-invalid-partial-specialization")
endif()
option(VQE_BUILD_TESTS "Build test programs" OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-invalid-partial-specialization")

if(APPLE)
  set(CMAKE_MACOSX_RPATH 1)
endif()

find_package(XACC REQUIRED)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${XACC_ROOT}" CACHE PATH "default install path" FORCE )
endif()

set(CppMicroServicesConfig_DIR "${XACC_ROOT}/share/cppmicroservices4/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${XACC_ROOT}/share/xacc")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CppMicroServicesConfig_DIR}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/tpls/backward-cpp")

include(CppMicroServicesConfig)

set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)

include_directories(${XACC_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/tpls/cppnumericalsolvers/cppoptlib)
include_directories(${CMAKE_SOURCE_DIR}/tpls/exprtk)

include_directories(${CMAKE_SOURCE_DIR}/tpls/backward-cpp)

link_directories(${XACC_LIBRARY_DIR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/mpi)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/ir)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/compiler)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/transformations)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/utils)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/task)

add_subdirectory(mpi)
add_subdirectory(ir)
add_subdirectory(transformations)
add_subdirectory(compiler)
add_subdirectory(task)
add_subdirectory(decorators)

if(PYTHON_INCLUDE_DIR)
   include_directories(${PYTHON_INCLUDE_DIR})
   add_subdirectory(python)
endif()

#add_executable(xacc-vqe xacc-vqe.cpp tpls/backward-cpp/backward.cpp)
#target_link_libraries(xacc-vqe xacc-vqe-ir xacc-vqe-tasks ${XACC_LIBRARIES} pthread) # petsc)
#if(APPLE)
#   set_target_properties(xacc-vqe PROPERTIES INSTALL_RPATH "@loader_path/../lib;@loader_path/../plugins")
#else()
#   set_target_properties(xacc-vqe PROPERTIES INSTALL_RPATH "$ORIGIN/../lib:$ORIGIN/../plugins")
#endif()

#install(TARGETS xacc-vqe DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

file(GLOB BENCHMARK ${CMAKE_CURRENT_SOURCE_DIR}/benchmark/*.py)
install(FILES ${BENCHMARK} DESTINATION ${CMAKE_INSTALL_PREFIX}/py-plugins)

file(GLOB PYDECORATORS ${CMAKE_CURRENT_SOURCE_DIR}/python/decorators/*.py)
install(FILES ${PYDECORATORS} DESTINATION ${CMAKE_INSTALL_PREFIX}/py-plugins)

configure_file("${CMAKE_SOURCE_DIR}/cmake/xaccvqe-config.cmake.in" "${CMAKE_BINARY_DIR}/xaccvqe-config.cmake" @ONLY)
install(FILES "${CMAKE_BINARY_DIR}/xaccvqe-config.cmake" DESTINATION ${CMAKE_INSTALL_PREFIX})
install( DIRECTORY "${CMAKE_SOURCE_DIR}/tpls/cppnumericalsolvers/cppoptlib" DESTINATION ${CMAKE_INSTALL_PREFIX}/include )
