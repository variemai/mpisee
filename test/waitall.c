#include <stdio.h>
#include <mpi.h>
#define NPROCS 4

int main (int argc, char *argv[]){
    int size;
    int rank;
    MPI_Comm newcomm;
    MPI_Request requests[2];
    MPI_Status status;
    int p;
    const int d = 10;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if ( size != NPROCS ){
        fprintf(stderr, "This app runs only with 4 ranks\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if ( rank % 2 == 0 ){
        MPI_Isend(&d, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, &requests[0]);
    }
    else
        MPI_Irecv(&p, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &requests[0]);


    MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank, &newcomm);
    MPI_Comm_rank(newcomm, &rank);
    if ( rank == 0 ){
        MPI_Isend(&rank, 1, MPI_INT, rank + 1, 0, newcomm, &requests[1]);
    }
    else
        MPI_Irecv(&p, 1, MPI_INT, rank - 1, 0, newcomm, &requests[1]);
    MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
    MPI_Finalize();

}
