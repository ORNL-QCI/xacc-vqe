project(xacc-vqe CXX)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_DISABLE_IN_SOURCE_BUILDS ON)
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")
include (tests)

find_package(PythonLibs REQUIRED)
find_package(Boost COMPONENTS filesystem program_options system python REQUIRED)
find_package(XACC REQUIRED)

include_directories(${PYTHON_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${XACC_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/tpls/cppnumericalsolvers/cppoptlib)
include_directories(${CMAKE_SOURCE_DIR}/vqe)
include_directories(${CMAKE_SOURCE_DIR}/vqe/ir)
include_directories(${CMAKE_SOURCE_DIR}/vqe/transformation/uccsd)
include_directories(${CMAKE_SOURCE_DIR}/vqe/compiler)
include_directories(${CMAKE_SOURCE_DIR}/tpls/eigen)

link_directories(${XACC_LIBRARY_DIR})

add_subdirectory(vqe)

add_executable(xacc-vqe main.cpp)
target_link_libraries(xacc-vqe ${SYMENGINE_LIBRARIES} ${Boost_LIBRARIES} xacc-vqe-ir ${XACC_LIBRARIES} pthread dl)
