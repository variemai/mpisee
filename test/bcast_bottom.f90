program mpi_bcast_bottom_example
    use mpi
    implicit none

    integer :: ierr, rank, size
    integer :: buffer(1)          ! Buffer for data to broadcast
    integer(kind=MPI_ADDRESS_KIND) :: addr    ! Base address of the buffer
    integer :: ones(1) = [1]      ! Block length for derived datatype
    integer :: bcast_type         ! Derived datatype

    ! Initialize MPI
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

    ! Get the base address of the buffer
    call MPI_Get_address(buffer, addr, ierr)

    ! Create a derived datatype for the buffer using MPI_BOTTOM
    call MPI_Type_create_struct(1, ones, [addr], [MPI_INTEGER], bcast_type, ierr)
    call MPI_Type_commit(bcast_type, ierr)

    ! Initialize the buffer differently on rank 0
    if (rank == 0) then
        buffer(1) = 42    ! Root initializes the buffer
    else
        buffer(1) = -1    ! Other ranks initialize the buffer to a dummy value
    end if

    ! Broadcast the data using MPI_BOTTOM
    call MPI_Bcast(MPI_BOTTOM, 1, bcast_type, 0, MPI_COMM_WORLD, ierr)

    ! Verify the result
    if (buffer(1) /= 42) then
        print *, "Error on rank ", rank, ": buffer(1) = ", buffer(1)
    else
        print *, "Rank ", rank, ": buffer(1) = ", buffer(1)
    end if

    ! Clean up
    call MPI_Type_free(bcast_type, ierr)
    call MPI_Finalize(ierr)
end program mpi_bcast_bottom_example
