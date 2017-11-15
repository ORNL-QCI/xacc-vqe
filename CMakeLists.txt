project(xacc-vqe CXX)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_DISABLE_IN_SOURCE_BUILDS ON)
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)

find_package(XACC REQUIRED)
set(CMAKE_INSTALL_RPATH "${XACC_ROOT}/lib;${XACC_ROOT}/lib/plugins/misc;/usr/local/lib64")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(CppMicroServicesConfig_DIR "${XACC_ROOT}/share/cppmicroservices3/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${XACC_ROOT}/share/xacc")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CppMicroServicesConfig_DIR}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

include(tests)
include(CppMicroServicesConfig)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(MPI REQUIRED)

find_package(Boost COMPONENTS filesystem program_options system mpi REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})

find_package(XACC REQUIRED)

include_directories(${PYTHON_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${XACC_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/tpls/cppnumericalsolvers/cppoptlib)
include_directories(${CMAKE_SOURCE_DIR}/tpls/exprtk)

link_directories(${XACC_LIBRARY_DIR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/ir)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/compiler)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/transformations)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/utils)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/task)


add_subdirectory(ir)
add_subdirectory(transformations)
add_subdirectory(compiler)
add_subdirectory(task)
add_subdirectory(python)

add_executable(xacc-vqe xacc-vqe.cpp)
target_link_libraries(xacc-vqe ${Boost_LIBRARIES} xacc-vqe-ir xacc-vqe-tasks ${XACC_LIBRARIES} pthread dl ${MPI_CXX_LIBRARIES}) # petsc)
