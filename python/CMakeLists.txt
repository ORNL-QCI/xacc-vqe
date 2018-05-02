include_directories(${XACC_ROOT}/include/pybind11/include)
include_directories(${PYTHON_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/task/tasks)

link_directories(${XACC_LIBRARY_DIR})

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-strict-aliasing -O2 -g -pipe -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -Wformat -fexceptions -fstack-protector-strong --param=ssp-buffer-size=4 -grecord-gcc-switches -m64 -mtune=generic -D_GNU_SOURCE -fPIC -fwrapv")

if(APPLE)
   set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
endif(APPLE)

add_library(pyxaccvqe SHARED xacc-vqe-py.cpp)
set_target_properties(pyxaccvqe PROPERTIES PREFIX "")
# do not mess with the ordering of libraries here!
target_link_libraries(pyxaccvqe xacc-vqe-tasks ${XACC_LIBRARIES})
if(APPLE)
   set_target_properties(pyxaccvqe PROPERTIES INSTALL_RPATH "@loader_path/lib;@loader_path/plugins")
   set_target_properties(pyxaccvqe PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
else()
   set_target_properties(pyxaccvqe PROPERTIES INSTALL_RPATH "$ORIGIN/lib:$ORIGIN/plugins")
   set_target_properties(pyxaccvqe PROPERTIES LINK_FLAGS "-shared")
endif()

install(TARGETS pyxaccvqe DESTINATION ${CMAKE_INSTALL_PREFIX})


file (GLOB PYTESTS tests/*.py)



