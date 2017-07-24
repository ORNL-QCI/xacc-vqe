set (PACKAGE_NAME "XACC VQE Runtime")
set (PACKAGE_DESCIPTION "XACC Variational Quantum Eigensolver Support")

set (LIBRARY_NAME xacc-vqe)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/ir)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/compiler)

file (GLOB_RECURSE HEADERS *.hpp)
file (GLOB SRC compiler/*.cpp ir/*.cpp)

add_library(${LIBRARY_NAME} SHARED ${SRC})

target_link_libraries(${LIBRARY_NAME} ${Boost_LIBRARIES} ${XACC_LIBRARIES} ${SYMENGINE_LIBRARIES})

add_subdirectory(ir)
add_subdirectory(compiler)

#install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION ${XACC_LIBRARY_DIR})