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

find_package(Boost COMPONENTS filesystem program_options system REQUIRED)
find_package(XACC REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${XACC_INCLUDE_DIRS})

add_subdirectory(vqe)
