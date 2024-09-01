cmake_minimum_required(VERSION 3.10)
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILER mpicxx)

add_executable(init init.cpp)
#add_executable(ctest0 ctest0.cpp)
#add_executable(ctest1 ctest1.cpp)
#add_executable(ctest_input ctest_input.cpp)
#add_executable(bcast_mt bcast_mt.cpp)
#add_executable(bcast_reduce_mt bcast_reduce_mt.cpp)
#add_executable(comm_create comm_create.cpp)
#add_executable(comm_split comm_split.cpp)
#add_executable(allg_rand_data allg_rand_data.cpp)
#add_executable(comm_free comm_free.cpp)

find_program(SRUN_EXEC srun)
if(SRUN_EXEC)
    message(STATUS "srun found at: ${SRUN_EXEC}")
    set(MPIEXEC ${SRUN_EXEC})
else()
    message(STATUS "srun not found in PATH")
endif()


add_test(NAME init COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${PROJECT_BINARY_DIR}/test/init)
add_test(NAME see_through_init COMMAND sh -c
    "python3 ${PROJECT_SOURCE_DIR}/mpisee-through/mpisee-through.py -i ${PROJECT_BINARY_DIR}/test/init.db | grep 'Init' > /dev/null 2>&1 && echo 'Init found' || (echo 'Init not found'; exit 1)" )
#add_test(NAME ctest0 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${PROJECT_BINARY_DIR}/test/ctest0)
#add_test(NAME ctest1 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${PROJECT_BINARY_DIR}/test/ctest1)
#add_test(NAME ctest_inputv COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${PROJECT_BINARY_DIR}/test/ctest_input -n 1032 -f file.txt -q abcdr -W -n 1032 -f file.txt -q abcdr -W -n 1032 -f file.txt -q abcdr -W)
#add_test(NAME bcast_mt COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${PROJECT_BINARY_DIR}/test/bcast_mt)
#add_test(NAME bcast_reduce_mt COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${PROJECT_BINARY_DIR}/test/bcast_reduce_mt)
#add_test(NAME comm_create COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${PROJECT_BINARY_DIR}/test/comm_create)
#add_test(NAME comm_split COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${PROJECT_BINARY_DIR}/test/comm_split)
#add_test(NAME allg_rand_data COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${PROJECT_BINARY_DIR}/test/allg_rand_data)
#add_test(NAME comm_free COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${PROJECT_BINARY_DIR}/test/comm_free)

set_tests_properties(init PROPERTIES ENVIRONMENT "MPISEE_OUTFILE=init.db")
set_tests_properties(see_through_init PROPERTIES DEPENDS init)

target_link_libraries(init PUBLIC mpisee)
#target_link_libraries(ctest0 PUBLIC mpisee)
#target_link_libraries(ctest1 PUBLIC mpisee)
#target_link_libraries(ctest_input PUBLIC mpisee)
#target_link_libraries(bcast_mt PUBLIC mpisee)
#target_link_libraries(bcast_reduce_mt PUBLIC mpisee)
#target_link_libraries(comm_create PUBLIC mpisee)
#target_link_libraries(comm_split PUBLIC mpisee)
#target_link_libraries(allg_rand_data PUBLIC mpisee)
#target_link_libraries(comm_free PUBLIC mpisee)

#set_property(TARGET ctest0 PROPERTY CROSSCOMPILING_EMULATOR '${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}')
# set(MPI_COMMAND mpirun -np 2 ctest0)
# add_custom_command(
#   TARGET ctest0
#   POST_BUILD
#   COMMAND ${MPI_CXXOMMAND})
#   set(MPI_CXXOMMAND mpirun -np 2 ctest1)
# add_custom_command(
#   TARGET ctest1
#   POST_BUILD
#   COMMAND ${MPI_CXXOMMAND})
#target_include_directories(init PRIVATE ${MPI_CXX_INCLUDE_PATH})
#target_compile_options(init PRIVATE ${MPI_CXX_COMPILE_FLAGS})
#target_link_libraries(init PRIVATE ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
#target_include_directories(init PUBLIC "${PROJECT_SOURCE_DIR}/")

#target_include_directories(ctest0 PRIVATE ${MPI_CXX_INCLUDE_PATH})
#target_compile_options(ctest0 PRIVATE ${MPI_CXX_COMPILE_FLAGS})
#target_link_libraries(ctest0 PRIVATE ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
#target_include_directories(ctest0 PUBLIC "${PROJECT_SOURCE_DIR}/")
#
#target_include_directories(ctest1 PRIVATE ${MPI_CXX_INCLUDE_PATH})
#target_compile_options(ctest1 PRIVATE ${MPI_CXX_COMPILE_FLAGS})
#target_link_libraries(ctest1 PRIVATE ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
#target_include_directories(ctest1 PUBLIC "${PROJECT_SOURCE_DIR}/")
#
#target_include_directories(ctest_input PRIVATE ${MPI_CXX_INCLUDE_PATH})
#target_compile_options(ctest_input PRIVATE ${MPI_CXX_COMPILE_FLAGS})
#target_link_libraries(ctest_input PRIVATE ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
#target_include_directories(ctest_input PUBLIC "${PROJECT_SOURCE_DIR}/")
#
#
#target_include_directories(bcast_mt PRIVATE ${MPI_CXX_INCLUDE_PATH})
#target_compile_options(bcast_mt PRIVATE ${MPI_CXX_COMPILE_FLAGS})
#target_link_libraries(bcast_mt PRIVATE ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
#target_include_directories(bcast_mt PUBLIC "${PROJECT_SOURCE_DIR}/")
#
#target_include_directories(bcast_reduce_mt PRIVATE ${MPI_CXX_INCLUDE_PATH})
#target_compile_options(bcast_reduce_mt PRIVATE ${MPI_CXX_COMPILE_FLAGS})
#target_link_libraries(bcast_reduce_mt PRIVATE ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
#target_include_directories(bcast_reduce_mt PUBLIC "${PROJECT_SOURCE_DIR}/")
#
#target_include_directories(comm_create PRIVATE ${MPI_CXX_INCLUDE_PATH})
#target_compile_options(comm_create PRIVATE ${MPI_CXX_COMPILE_FLAGS})
#target_link_libraries(comm_create PRIVATE ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
#target_include_directories(comm_create PUBLIC "${PROJECT_SOURCE_DIR}/")
#
#target_include_directories(comm_split PRIVATE ${MPI_CXX_INCLUDE_PATH})
#target_compile_options(comm_split PRIVATE ${MPI_CXX_COMPILE_FLAGS})
#target_link_libraries(comm_split PRIVATE ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
#target_include_directories(comm_split PUBLIC "${PROJECT_SOURCE_DIR}/")
#
#target_include_directories(allg_rand_data PRIVATE ${MPI_CXX_INCLUDE_PATH})
#target_compile_options(allg_rand_data PRIVATE ${MPI_CXX_COMPILE_FLAGS})
#target_link_libraries(allg_rand_data PRIVATE ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
#target_include_directories(allg_rand_data PUBLIC "${PROJECT_SOURCE_DIR}/")
#
#target_include_directories(comm_free PRIVATE ${MPI_CXX_INCLUDE_PATH})
#target_compile_options(comm_free PRIVATE ${MPI_CXX_COMPILE_FLAGS})
#target_link_libraries(comm_free PRIVATE ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
#target_include_directories(comm_free PUBLIC "${PROJECT_SOURCE_DIR}/")
