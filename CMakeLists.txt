cmake_minimum_required(VERSION 3.10)
project(mpisee VERSION 1.0 DESCRIPTION "MPI Communicator Profiler")
find_package(MPI REQUIRED)
message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPI_C_COMPILER_INCLUDE_DIRS}")

add_subdirectory(datastructlib)

add_library(mpisee SHARED
  commprof.c
  utils.c
  )
target_include_directories(mpisee PUBLIC "$(PROJECT_SOURCE_DIR)/datastructlib")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -I./datastructlib")
configure_file(commprofConfig.h.in commprofConfig.h)
target_include_directories(mpisee PRIVATE ${MPI_C_INCLUDE_PATH})
target_compile_options(mpisee PRIVATE ${MPI_C_COMPILE_FLAGS})
#target_link_libraries(mpisee PUBLIC ${MPI_C_LIB_NAMES})
target_link_libraries(mpisee PRIVATE ${MPI_C_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
target_link_libraries(mpisee PUBLIC datastructlib)
set_target_properties(mpisee PROPERTIES SOVERSION 1)
add_subdirectory(test)
#add_test(Barrier )


#target_include_directories(my_mpi_bin PRIVATE ${MPI_CXX_INCLUDE_PATH})
#target_compile_options(my_mpi_bin PRIVATE ${MPI_CXX_COMPILE_FLAGS} my_extra_compile_flags)
#target_link_libraries(my_mpi_bin ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_F
