program allreduce_bottom_example
    use mpi
    implicit none

    integer :: ierr, rank, size, i
    integer, dimension(:), allocatable :: my_values, result_values
    integer(kind=MPI_ADDRESS_KIND) :: disp_send, disp_recv
    integer(kind=MPI_ADDRESS_KIND), dimension(2) :: displacements
    integer, dimension(2) :: block_lengths
    integer :: datatype
    integer(kind=MPI_ADDRESS_KIND) :: base_address

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

    ! Each process has a unique array of values
    allocate(my_values(size))
    allocate(result_values(size))

    ! Initialize the arrays
    my_values = rank + 1
    result_values = 0

    ! Determine the base address of memory
    call MPI_Get_address(my_values, base_address, ierr)
    call MPI_Get_address(my_values(1), disp_send, ierr)
    call MPI_Get_address(result_values, disp_recv, ierr)

    ! Calculate displacements relative to MPI_BOTTOM
    displacements(1) = disp_send - base_address
    displacements(2) = disp_recv - base_address

    block_lengths = [size, size]  ! The block sizes for both buffers

    ! Create a derived datatype for the operation
    call MPI_Type_create_hindexed(2, block_lengths, displacements, MPI_INTEGER, datatype, ierr)
    call MPI_Type_commit(datatype, ierr)

    ! Perform the reduction operation using MPI_BOTTOM
    call MPI_Allreduce(MPI_BOTTOM, MPI_BOTTOM, 1, datatype, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Check the result
    print *, 'Process', rank, 'my_values:', my_values
    print *, 'Process', rank, 'result_values after MPI_Allreduce with MPI_BOTTOM:', result_values

    ! Clean up
    call MPI_Type_free(datatype, ierr)
    deallocate(my_values)
    deallocate(result_values)
    call MPI_Finalize(ierr)
end program allreduce_bottom_example
