cmake_minimum_required(VERSION 3.10)
find_package(MPI)
set(CMAKE_CXX_COMPILER mpicxx)

add_executable(init.out init.cpp)
add_executable(allgatherv_in_place.out allgatherv_inplace.f90)
add_executable(allgather_in_place.out allgather_inplace.f90)
add_executable(allreduce_in_place.out allreduce_inplace.f90)
add_executable(reduce_in_place.out reduce_inplace.f90)
add_executable(bcast_bottom.out bcast_bottom.f90)
add_executable(alltoallw_simple.out alltoallw_simple.f90)


#find_program(SRUN_EXEC srun)
#if(SRUN_EXEC)
#    message(STATUS "srun found at: ${SRUN_EXEC}")
#    set(MPIEXEC ${SRUN_EXEC})
#else()
#    message(STATUS "srun not found in PATH")
#endif()


add_test(NAME init COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${PROJECT_BINARY_DIR}/test/init.out)
add_test(NAME see_through_init COMMAND sh -c
    "python3 ${PROJECT_SOURCE_DIR}/mpisee-through/mpisee-through.py -i ${PROJECT_BINARY_DIR}/test/init.db | grep 'Init' > /dev/null 2>&1 && echo 'Init found' || (echo 'Init not found'; exit 1)" )
add_test(NAME allgather_in_place COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${PROJECT_BINARY_DIR}/test/allgather_in_place.out)
add_test(NAME see_through_allgather_in_place COMMAND sh -c
    "python3 ${PROJECT_SOURCE_DIR}/mpisee-through/mpisee-through.py -i ${PROJECT_BINARY_DIR}/test/allgather_in_place.db | grep 'Allgather' > /dev/null 2>&1 && echo 'Allgather found' || (echo 'Allgather not found'; exit 1)" )
add_test(NAME allgatherv_in_place COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${PROJECT_BINARY_DIR}/test/allgatherv_in_place.out)
add_test(NAME see_through_allgatherv_in_place COMMAND sh -c
    "python3 ${PROJECT_SOURCE_DIR}/mpisee-through/mpisee-through.py -i ${PROJECT_BINARY_DIR}/test/allgatherv_in_place.db | grep 'Allgatherv' > /dev/null 2>&1 && echo 'Allgatherv found' || (echo 'Allgatherv not found'; exit 1)" )
add_test(NAME allreduce_in_place COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${PROJECT_BINARY_DIR}/test/allreduce_in_place.out)
add_test(NAME see_through_allreduce_in_place COMMAND sh -c
    "python3 ${PROJECT_SOURCE_DIR}/mpisee-through/mpisee-through.py -i ${PROJECT_BINARY_DIR}/test/allreduce_in_place.db | grep 'Allreduce' > /dev/null 2>&1 && echo 'Allreduce found' || (echo 'Allreduce not found'; exit 1)" )
add_test(NAME reduce_in_place COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${PROJECT_BINARY_DIR}/test/reduce_in_place.out)
add_test(NAME see_through_reduce_in_place COMMAND sh -c
    "python3 ${PROJECT_SOURCE_DIR}/mpisee-through/mpisee-through.py -i ${PROJECT_BINARY_DIR}/test/reduce_in_place.db | grep 'Reduce' > /dev/null 2>&1 && echo 'Reduce found' || (echo 'Reduce not found'; exit 1)" )
add_test(NAME bcast_bottom COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${PROJECT_BINARY_DIR}/test/bcast_bottom.out)
add_test(NAME see_through_bcast_bottom COMMAND sh -c "python3 ${PROJECT_SOURCE_DIR}/mpisee-through/mpisee-through.py -i ${PROJECT_BINARY_DIR}/test/bcast_bottom.db | grep 'Bcast' > /dev/null 2>&1 && echo 'Bcast found' || (echo 'Bcast not found'; exit 1)" )
add_test(NAME alltoallw_simple COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${PROJECT_BINARY_DIR}/test/alltoallw_simple.out)
add_test(NAME see_through_alltoallw_simple COMMAND sh -c "python3 ${PROJECT_SOURCE_DIR}/mpisee-through/mpisee-through.py -i ${PROJECT_BINARY_DIR}/test/alltoallw_simple.db | grep 'Alltoallw' > /dev/null 2>&1 && echo 'Alltoallw found' || (echo 'Alltoallw not found'; exit 1)" )

# Always run cleanup as the last test
add_test(NAME cleanup
        COMMAND sh -c "echo 'Cleaning up test artifacts...' && rm -f ${PROJECT_BINARY_DIR}/test/*.db"
)


set_tests_properties(init PROPERTIES ENVIRONMENT "MPISEE_OUTFILE=init.db")
set_tests_properties(allgather_in_place PROPERTIES ENVIRONMENT "MPISEE_OUTFILE=allgather_in_place.db")
set_tests_properties(allgatherv_in_place PROPERTIES ENVIRONMENT "MPISEE_OUTFILE=allgatherv_in_place.db")
set_tests_properties(allreduce_in_place PROPERTIES ENVIRONMENT "MPISEE_OUTFILE=allreduce_in_place.db")
set_tests_properties(reduce_in_place PROPERTIES ENVIRONMENT "MPISEE_OUTFILE=reduce_in_place.db")
set_tests_properties(bcast_bottom PROPERTIES ENVIRONMENT "MPISEE_OUTFILE=bcast_bottom.db")
set_tests_properties(alltoallw_simple PROPERTIES ENVIRONMENT "MPISEE_OUTFILE=alltoallw_simple.db")

set_tests_properties(see_through_init PROPERTIES DEPENDS init)
set_tests_properties(see_through_allgather_in_place PROPERTIES DEPENDS allgather_in_place)
set_tests_properties(see_through_allgatherv_in_place PROPERTIES DEPENDS allgatherv_in_place)
set_tests_properties(see_through_allreduce_in_place PROPERTIES DEPENDS allreduce_in_place)
set_tests_properties(see_through_reduce_in_place PROPERTIES DEPENDS reduce_in_place)
set_tests_properties(see_through_bcast_bottom PROPERTIES DEPENDS bcast_bottom)
set_tests_properties(see_through_alltoallw_simple PROPERTIES DEPENDS alltoallw_simple)


target_link_libraries(init.out PUBLIC mpisee)
target_link_libraries(allgather_in_place.out PUBLIC mpisee)
target_link_libraries(allgatherv_in_place.out PUBLIC mpisee)
target_link_libraries(allreduce_in_place.out PUBLIC mpisee)
target_link_libraries(reduce_in_place.out PUBLIC mpisee)
target_link_libraries(bcast_bottom.out PUBLIC mpisee)
target_link_libraries(alltoallw_simple.out PUBLIC mpisee)

add_custom_target(manual_cleanup
        COMMAND ${CMAKE_COMMAND} -E rm -f ${PROJECT_BINARY_DIR}/test/*.db
        COMMAND ${CMAKE_COMMAND} -E rm -f ${PROJECT_BINARY_DIR}/test/*.out
        # Add more paths or file patterns as needed
        COMMENT "Cleaning up test artifacts..."
)

# Combined target to run tests and then clean up
add_custom_target(run_and_cleanup_tests
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target cleanup
        COMMENT "Running tests and cleaning up..."
)


target_include_directories(init.out PRIVATE ${MPI_CXX_INCLUDE_PATH})
#target_compile_options(init PRIVATE ${MPI_CXX_COMPILE_FLAGS})
target_link_libraries(init.out PRIVATE ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
#target_include_directories(init PUBLIC "${PROJECT_SOURCE_DIR}/")

target_include_directories(allreduce_in_place.out PRIVATE ${MPI_Fortran_INCLUDE_PATH})
#target_compile_options(init PRIVATE ${MPI_CXX_COMPILE_FLAGS})
target_link_libraries(allreduce_in_place.out PRIVATE ${MPI_Fortran_LIBRARIES} ${MPI_Fortran_LINK_FLAGS})

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
