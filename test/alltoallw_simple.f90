PROGRAM main
    USE mpi
    IMPLICIT NONE

    INTEGER :: size
    INTEGER :: rank
    INTEGER :: ierr
    INTEGER :: i
    INTEGER :: int_size

    INTEGER, ALLOCATABLE :: recvbuf(:)
    INTEGER, ALLOCATABLE :: recvcounts(:)
    INTEGER, ALLOCATABLE :: recvdispls(:)
    INTEGER, ALLOCATABLE :: recvtypes(:)

    INTEGER, ALLOCATABLE :: sendbuf(:)
    INTEGER, ALLOCATABLE :: sendcounts(:)
    INTEGER, ALLOCATABLE :: senddispls(:)
    INTEGER, ALLOCATABLE :: sendtypes(:)

    CALL MPI_Init(ierr)
    CALL MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
    CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    CALL MPI_Type_size(MPI_INTEGER, int_size, ierr)

    ALLOCATE(recvbuf(size))
    ALLOCATE(recvcounts(size))
    ALLOCATE(recvdispls(size))
    ALLOCATE(recvtypes(size))

    ALLOCATE(sendbuf(size))
    ALLOCATE(sendcounts(size))
    ALLOCATE(senddispls(size))
    ALLOCATE(sendtypes(size))
    DO i = 1, size
        recvcounts(i) = 1
        recvdispls(i) = (i-1) * int_size
        recvtypes(i)  = MPI_INTEGER
        recvbuf(i)    = rank
        sendcounts(i) = 1
        senddispls(i) = (i-1) * int_size
        sendtypes(i)  = MPI_INTEGER
        sendbuf(i)    = rank*10+1
    END DO

    CALL MPI_Alltoallw(sendbuf, sendcounts, senddispls, sendtypes, &
                       recvbuf, recvcounts, recvdispls, recvtypes, &
                       MPI_COMM_WORLD, ierr)

    ! Print receive buffer
    PRINT *, "RANK:", rank, "Received buffer:", recvbuf

    CALL MPI_Finalize(ierr)
END PROGRAM
