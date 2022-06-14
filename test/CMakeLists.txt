cmake_minimum_required(VERSION 3.10)
find_package(MPI REQUIRED)

add_executable(ctest0 ctest0.c)
add_executable(ctest1 ctest1.c)
add_executable(ctest_input ctest_input.c)
add_executable(ctest_inputv ctest_inputv.c)
add_executable(bcast_mt bcast_mt.c)
add_executable(bcast_reduce_mt bcast_reduce_mt.c)
add_executable(comm_create comm_create.c)
add_executable(comm_split comm_split.c)
add_executable(allg_rand_data allg_rand_data.c)
add_executable(comm_free comm_free.c)


add_test(NAME ctest0 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${PROJECT_BINARY_DIR}/test/ctest0)
add_test(NAME ctest1 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${PROJECT_BINARY_DIR}/test/ctest1)
add_test(NAME ctest_input COMMAND  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${PROJECT_BINARY_DIR}/test/ctest_input -n 1032 -f file.txt -q abcdr -W)
add_test(NAME ctest_inputv COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${PROJECT_BINARY_DIR}/test/ctest_inputv)
add_test(NAME bcast_mt COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${PROJECT_BINARY_DIR}/test/bcast_mt)
add_test(NAME bcast_reduce_mt COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${PROJECT_BINARY_DIR}/test/bcast_reduce_mt)
add_test(NAME comm_create COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${PROJECT_BINARY_DIR}/test/comm_create)
add_test(NAME comm_split COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${PROJECT_BINARY_DIR}/test/comm_split)
add_test(NAME allg_rand_data COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${PROJECT_BINARY_DIR}/test/allg_rand_data)
add_test(NAME comm_free COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${PROJECT_BINARY_DIR}/test/comm_free)

target_link_libraries(ctest0 PUBLIC mpisee)
target_link_libraries(ctest1 PUBLIC mpisee)
target_link_libraries(ctest_input PUBLIC mpisee)
target_link_libraries(ctest_inputv PUBLIC mpisee)
target_link_libraries(bcast_mt PUBLIC mpisee)
target_link_libraries(bcast_reduce_mt PUBLIC mpisee)
target_link_libraries(comm_create PUBLIC mpisee)
target_link_libraries(comm_split PUBLIC mpisee)
target_link_libraries(allg_rand_data PUBLIC mpisee)
target_link_libraries(comm_free PUBLIC mpisee)

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

target_include_directories(ctest_inputv PRIVATE ${MPI_C_INCLUDE_PATH})
target_compile_options(ctest_inputv PRIVATE ${MPI_C_COMPILE_FLAGS})
target_link_libraries(ctest_inputv PRIVATE ${MPI_C_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
target_include_directories(ctest_inputv PUBLIC "${PROJECT_SOURCE_DIR}/")

target_include_directories(bcast_mt PRIVATE ${MPI_C_INCLUDE_PATH})
target_compile_options(bcast_mt PRIVATE ${MPI_C_COMPILE_FLAGS})
target_link_libraries(bcast_mt PRIVATE ${MPI_C_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
target_include_directories(bcast_mt PUBLIC "${PROJECT_SOURCE_DIR}/")

target_include_directories(bcast_reduce_mt PRIVATE ${MPI_C_INCLUDE_PATH})
target_compile_options(bcast_reduce_mt PRIVATE ${MPI_C_COMPILE_FLAGS})
target_link_libraries(bcast_reduce_mt PRIVATE ${MPI_C_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
target_include_directories(bcast_reduce_mt PUBLIC "${PROJECT_SOURCE_DIR}/")

target_include_directories(comm_create PRIVATE ${MPI_C_INCLUDE_PATH})
target_compile_options(comm_create PRIVATE ${MPI_C_COMPILE_FLAGS})
target_link_libraries(comm_create PRIVATE ${MPI_C_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
target_include_directories(comm_create PUBLIC "${PROJECT_SOURCE_DIR}/")

target_include_directories(comm_split PRIVATE ${MPI_C_INCLUDE_PATH})
target_compile_options(comm_split PRIVATE ${MPI_C_COMPILE_FLAGS})
target_link_libraries(comm_split PRIVATE ${MPI_C_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
target_include_directories(comm_split PUBLIC "${PROJECT_SOURCE_DIR}/")

target_include_directories(allg_rand_data PRIVATE ${MPI_C_INCLUDE_PATH})
target_compile_options(allg_rand_data PRIVATE ${MPI_C_COMPILE_FLAGS})
target_link_libraries(allg_rand_data PRIVATE ${MPI_C_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
target_include_directories(allg_rand_data PUBLIC "${PROJECT_SOURCE_DIR}/")

target_include_directories(comm_free PRIVATE ${MPI_C_INCLUDE_PATH})
target_compile_options(comm_free PRIVATE ${MPI_C_COMPILE_FLAGS})
target_link_libraries(comm_free PRIVATE ${MPI_C_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
target_include_directories(comm_free PUBLIC "${PROJECT_SOURCE_DIR}/")
