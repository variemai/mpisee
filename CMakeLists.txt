cmake_minimum_required(VERSION 3.10)
project(mpisee VERSION 4.0 DESCRIPTION "mpisee: A communicator-centric profiler" LANGUAGES CXX Fortran)

find_package(MPI)
if(NOT MPI_FOUND)
    message(WARNING "MPI not found. mpisee will not compile.")
endif()
find_package(SQLite3)
if(NOT SQLite3_FOUND)
    message(WARNING "SQLite3 not found. mpisee will not compile.")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -Wall")

# Set default values
set(NUM_BUCKETS 8 CACHE STRING "Number of buckets")
set(BUCKETS "128,1024,8192,65536,262144,1048576,33554432" CACHE STRING "Message size buckets in bytes")

# Validation Checks
if(NOT NUM_BUCKETS MATCHES "^[1-9][0-9]*$")
    message(FATAL_ERROR "NUM_BUCKETS must be a positive integer")
endif()
if(NUM_BUCKETS LESS 2)
    message(FATAL_ERROR "NUM_BUCKETS must be at least 2 for meaningful bucketing")
endif()
string(REGEX MATCHALL "[0-9]+" bucket_values "${BUCKETS}")
list(LENGTH bucket_values num_values)
math(EXPR NB "${NUM_BUCKETS} - 1")
if(NOT num_values EQUAL NB)
    message(FATAL_ERROR "The number of entries in BUCKETS must be one less than NUM_BUCKETS")
endif()

# Extract major and minor version numbers
string(REGEX MATCH "^([0-9]+)\\.([0-9]+)" _ ${PROJECT_VERSION})
set(MPI_SEE_MAJOR_VERSION ${CMAKE_MATCH_1})
set(MPI_SEE_MINOR_VERSION ${CMAKE_MATCH_2})

configure_file(${CMAKE_SOURCE_DIR}/utils.h.in ${CMAKE_SOURCE_DIR}/include/utils.h)
message(STATUS "Generated utils.h in ${CMAKE_SOURCE_DIR}/include")

add_library(mpisee SHARED
        src/commprof.cpp
        src/pt2pt_wrappers.cpp
        src/collectives_wrappers.cpp
        src/neighbor_wrappers.cpp
        src/onesided_wrappers.cpp
        src/utils.cpp
        src/create_db.cpp
        src/mpisee_fortran.cpp
        src/mpisee_fortran_init.f90
)

target_include_directories(mpisee PRIVATE ${MPI_CXX_INCLUDE_DIRS} ${SQLite3_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR})
target_link_libraries(mpisee PRIVATE ${MPI_CXX_LIBRARIES} ${SQLite3_LIBRARIES})

set_target_properties(mpisee PROPERTIES SOVERSION 1)

enable_testing()
add_subdirectory(test)
find_program(SRUN_EXEC srun)

if(SRUN_EXEC)
    set(MPIEXEC "${SRUN_EXEC}")
    message(STATUS "Found srun: ${SRUN_EXEC}")
else()
    message(STATUS "srun not found, trying MPIEXEC")
    find_program(MPIEXEC mpirun)
    if (MPIEXEC)
        message(STATUS "MPIEXEC found at: ${MPIEXEC}")
    else()
        message(FATAL_ERROR "MPIEXEC not found. Please install MPI and set MPIEXEC to the MPI launcher")
    endif ()
endif()

# Installation
set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Installation prefix")
install(TARGETS mpisee
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin)

# Maybe in the future
# install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/ DESTINATION include)
install(FILES "${CMAKE_SOURCE_DIR}/mpisee-through/mpisee-through.py" DESTINATION bin
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE)
