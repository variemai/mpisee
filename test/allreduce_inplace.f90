program allreduce_in_place_example
    use mpi
    implicit none

    integer :: ierr, rank, size
    integer :: my_value, result_value

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

    ! Initialize each process with a unique value
    my_value = rank + 1

    ! Use MPI_IN_PLACE to perform in-place reduction
    result_value = my_value
    call MPI_Allreduce(MPI_IN_PLACE, result_value, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Each process should now have the sum of all values
    print *, 'Process', rank, 'result after MPI_Allreduce with MPI_IN_PLACE:', result_value

    call MPI_Finalize(ierr)
end program allreduce_in_place_example
