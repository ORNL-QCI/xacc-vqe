find_package(SLEPc) 

if (SLEPC_FOUND) 
   add_subdirectory(slepc)
endif()

add_subdirectory(tao-vqe)
