cmake_minimum_required(VERSION 3.10)
project(mpisee VERSION 3.3 DESCRIPTION "mpisee: A communicator-centric profiler")

find_package(MPI)
if(NOT MPI_FOUND)
   message(WARNING "MPI not found. mpisee will not compile.")
endif()
find_package(SQLite3)
if(NOT SQLite3_FOUND)
   message(WARNING "SQLite3 not found. mpisee will not compile.")
endif()


set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -Wall")
# Set default values
set(NUM_BUCKETS 8 CACHE STRING "Number of buckets")
set(BUCKETS "128,1024,8192,65536,262144,1048576,33554432" CACHE STRING "Message size buckets in bytes")

# Validation Checks

# 1. Check if NUM_BUCKETS is a positive integer
if(NOT NUM_BUCKETS MATCHES "^[1-9][0-9]*$")
  message(FATAL_ERROR "NUM_BUCKETS must be a positive integer")
endif()
# 2. Check if NUM_BUCKETS is at least 2
if(NUM_BUCKETS LESS 2)
    message(FATAL_ERROR "NUM_BUCKETS must be at least 2 for meaningful bucketing")
endif()
# 2. Check if the specified BUCKETS values are integers:
string(REGEX MATCHALL "[0-9]+" bucket_values "${BUCKETS}")
list(LENGTH bucket_values num_values)
math(EXPR NB "${NUM_BUCKETS} - 1")  # Corrected calculation
# 3. Check if the number of entries in BUCKETS matches NUM_BUCKETS-1
if(NOT num_values EQUAL NB )
    message(STATUS "Lengh of bucket_values: ${num_values}, NUM_BUCKETS: ${NUM_BUCKETS}")
    message(FATAL_ERROR "The number of entries in BUCKETS must be one less than NUM_BUCKETS")
endif()

# Extract major and minor version numbers
string(REGEX MATCH "^([0-9]+)\\.([0-9]+)" _ ${PROJECT_VERSION})
set(MPI_SEE_MAJOR_VERSION ${CMAKE_MATCH_1})
set(MPI_SEE_MINOR_VERSION ${CMAKE_MATCH_2})
# Generate the header file (the configure_file part remains the same)
configure_file(utils.h.in ${CMAKE_SOURCE_DIR}/utils.h)
message(STATUS "Generated utils.h in ${CMAKE_SOURCE_DIR}")


add_library(mpisee SHARED commprof.cpp pt2pt_wrappers.cpp collectives_wrappers.cpp neighbor_wrappers.cpp onesided_wrappers.cpp utils.cpp create_db.cpp)
#target_include_directories(create_db PRIVATE ${MPI_CXX_INCLUDE_DIRS})
#target_link_libraries(create_db PRIVATE ${MPI_CXX_LIBRARIES})
#target_include_directories(mpisee PRIVATE ${MPI_CXX_INCLUDE_DIRS})
#target_link_libraries(mpisee PRIVATE ${MPI_CXX_LIBRARIES})
target_link_libraries(mpisee PRIVATE ${SQLite3_LIBRARIES})
target_include_directories(mpisee PRIVATE ${SQLite3_INCLUDE_DIRS})

#if(MPI_CXX_COMPILE_FLAGS)
#    target_compile_options(mpisee PRIVATE ${MPI_CXX_COMPILE_FLAGS})
#endif()
#
#if(MPI_CXX_LINK_FLAGS)
#    set_target_properties(mpisee PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
#endif()

set_target_properties(mpisee PROPERTIES SOVERSION 1)
#message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPI_CXX_INCLUDE_DIRS}")
message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

enable_testing()
add_subdirectory(test)
find_program(SRUN_EXEC srun)

if(SRUN_EXEC)
    set(MPIEXEC "${SRUN_EXEC}")
    message(STATUS "Found srun: ${SRUN_EXEC}")
else()
    message(STATUS "srun not found, using default MPIEXEC")
endif()

# Installation
set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Installation prefix")
install(TARGETS mpisee
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin)

#install(FILES "${CMAKE_SOURCE_DIR}/utils.h" DESTINATION include)
install(FILES "${CMAKE_SOURCE_DIR}/mpisee-through/mpisee-through.py" DESTINATION bin
          PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                  GROUP_READ GROUP_EXECUTE
                  WORLD_READ WORLD_EXECUTE)
