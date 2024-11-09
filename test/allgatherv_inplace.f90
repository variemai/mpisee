program allgatherv_inplace_test
    use mpi
    implicit none

    integer :: ierr, rank, size
    integer, parameter :: n = 2  ! Assuming 2 processes
    integer :: recvbuf(n), recvcounts(n), displs(n)
    integer :: i

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

    ! Initialize each rank's data in recvbuf (used in place of sendbuf)
    recvbuf = 0  ! Initialize all elements to zero
    recvbuf(rank + 1) = rank + 1

    ! Initialize recvcounts and displacements for MPI_Allgatherv
    do i = 1, n
        recvcounts(i) = 1
        displs(i) = i - 1
    end do

    if (rank == 0) then
        print *, 'Before MPI_Allgatherv with MPI_IN_PLACE, recvbuf =', recvbuf
    end if

    ! Using MPI_IN_PLACE for MPI_Allgatherv
    call MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                        recvbuf, recvcounts, displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    print *, 'Rank ', rank, ', recvbuf after MPI_Allgatherv with MPI_IN_PLACE: ', recvbuf

    call MPI_Finalize(ierr)
end program allgatherv_inplace_test
