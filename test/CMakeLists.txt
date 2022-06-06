cmake_minimum_required(VERSION 3.10)
find_package(MPI REQUIRED)

add_executable(ctest0 ctest0.c)
add_executable(ctest1 ctest1.c)
add_executable(ctest_input ctest_input.c)
add_test(NAME ctest0 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${PROJECT_BINARY_DIR}/test/ctest0)
add_test(NAME ctest1 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${PROJECT_BINARY_DIR}/test/ctest1)
add_test(NAME ctest_input COMMAND ${PROJECT_BINARY_DIR}/test/ctest_input -n 1032 -f file.txt -q abcdr -W)
target_link_libraries(ctest0 PUBLIC mpisee)
target_link_libraries(ctest1 PUBLIC mpisee)
target_link_libraries(ctest_input PUBLIC mpisee)
message(STATUS "Test: ${PROJECT_BINARY_DIR}")
#set_property(TARGET ctest0 PROPERTY CROSSCOMPILING_EMULATOR '${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}')
# set(MPI_COMMAND mpirun -np 2 ctest0)
# add_custom_command(
#   TARGET ctest0
#   POST_BUILD
#   COMMAND ${MPI_COMMAND})
#   set(MPI_COMMAND mpirun -np 2 ctest1)
# add_custom_command(
#   TARGET ctest1
#   POST_BUILD
#   COMMAND ${MPI_COMMAND})
target_include_directories(ctest0 PRIVATE ${MPI_C_INCLUDE_PATH})
target_compile_options(ctest0 PRIVATE ${MPI_C_COMPILE_FLAGS})
target_link_libraries(ctest0 PRIVATE ${MPI_C_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
target_include_directories(ctest0 PUBLIC "${PROJECT_SOURCE_DIR}/")

target_include_directories(ctest1 PRIVATE ${MPI_C_INCLUDE_PATH})
target_compile_options(ctest1 PRIVATE ${MPI_C_COMPILE_FLAGS})
target_link_libraries(ctest1 PRIVATE ${MPI_C_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
target_include_directories(ctest1 PUBLIC "${PROJECT_SOURCE_DIR}/")

target_include_directories(ctest_input PRIVATE ${MPI_C_INCLUDE_PATH})
target_compile_options(ctest_input PRIVATE ${MPI_C_COMPILE_FLAGS})
target_link_libraries(ctest_input PRIVATE ${MPI_C_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
target_include_directories(ctest_input PUBLIC "${PROJECT_SOURCE_DIR}/")
