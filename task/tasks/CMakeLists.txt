#find_package(Boost COMPONENTS mpi)
#find_package(MPI)
#find_package(PETSc) 

if (Boost_FOUND AND MPI_CXX_FOUND AND PETSC_FOUND) 
   add_subdirectory(petsc)
endif()
