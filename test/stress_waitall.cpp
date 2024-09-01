#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#define ITERS 512              /* Default number of iterations */

int main (int argc, char *argv[]){
    int size, rank,iters,i,color,t;
    MPI_Comm *comms;
    MPI_Request *requests;
    MPI_Status status;
    MPI_Init(&argc,&argv);
    if ( argc < 2 ){
        iters = ITERS;
    }
    else{
        iters = atoi(argv[1]);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if ( rank % 2 == 0 )
        color = rank;
    else
        color = rank - 1;
    comms = (MPI_Comm*) malloc ( sizeof(MPI_Comm)*iters );
    requests = (MPI_Request*) malloc ( sizeof(MPI_Request)*iters );
    for ( i = 0; i<iters; i++ ){
        MPI_Comm_split(MPI_COMM_WORLD, color, rank % 2, &comms[i]);
        if ( rank % 2 == 0 )
            MPI_Isend(&rank,  1, MPI_INT, 1, 0, comms[i], &requests[i]);
        else
            MPI_Irecv(&t, 1, MPI_INT, 0, 0, comms[i], &requests[i]);
    }
    MPI_Waitall(iters,requests, MPI_STATUSES_IGNORE);
    for ( i = 0; i<iters; i++ ){
        MPI_Comm_free(&comms[i]);
    }
    MPI_Finalize();
    return 0;
}
