program reduce_inplace_test
    use mpi
    implicit none
    integer :: ierr, rank, size
    integer :: recvbuf

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

    ! Initialize recvbuf with each rank's data
    recvbuf = rank + 1

    if (rank == 0) then
        print *, "Before MPI_Reduce with MPI_IN_PLACE, recvbuf = ", recvbuf
    end if

    ! Root process uses MPI_IN_PLACE, others use recvbuf as sendbuf
    if (rank == 0) then
        call MPI_Reduce(MPI_IN_PLACE, recvbuf, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    else
        call MPI_Reduce(recvbuf, recvbuf, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    endif

    if (rank == 0) then
        print *, "After MPI_Reduce with MPI_IN_PLACE, recvbuf = ", recvbuf
    end if

    call MPI_Finalize(ierr)
end program reduce_inplace_test
