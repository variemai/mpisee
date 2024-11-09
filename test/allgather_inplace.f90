program allgather_test
    use mpi
    implicit none

    integer :: ierr, rank, size
    integer, parameter :: n = 2  ! Assuming 2 processes for simplicity
    integer :: recvbuf(n)
    integer :: i

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

    ! Initialize each rank's data
    recvbuf = rank + 1

    if (rank == 0) then
        print *, 'Before MPI_Allgather, recvbuf =', recvbuf
    end if

    ! Perform MPI_Allgather with MPI_IN_PLACE and MPI_DATATYPE_NULL for send type
    call MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                       recvbuf, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    ! Display the gathered data
    print *, 'Rank ', rank, ', recvbuf after MPI_Allgather: ', recvbuf

    call MPI_Finalize(ierr)
end program allgather_test
