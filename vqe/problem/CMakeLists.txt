set (PACKAGE_NAME "XACC VQE Problem")
set (PACKAGE_DESCIPTION "")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Gather tests
file (GLOB test_files tests/*.cpp)
add_tests("${test_files}" "${CMAKE_CURRENT_SOURCE_DIR}" "${Boost_LIBRARIES};pthread;dl;xacc-vqe-ir;xacc-vqe-irtransformations")

