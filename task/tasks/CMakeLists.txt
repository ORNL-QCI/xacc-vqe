find_package(Boost COMPONENTS mpi)
find_package(MPI)

if (BOOST_FOUND AND MPI_FOUND) 
   add_subdirectory(slepc)
endif()
