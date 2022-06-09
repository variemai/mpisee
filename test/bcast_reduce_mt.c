#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "../utils.h"
#include "../commprof.h"
#include <assert.h>

/* Metamorphic relations between Bcat and Reduce functions */

int main (int argc, char *argv[]){
    int flag,rank;
    MPI_Comm splitcomm;
    int buffer;
    prof_attrs *communicator = NULL;
    flag = 0;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    MPI_Comm_split(MPI_COMM_WORLD,rank % 2,rank,&splitcomm);
    /* 1st Bcast call with 2 MPI_INT */
    MPI_Bcast(&rank, 1, MPI_INT, 0, splitcomm);
    MPI_Reduce(&rank, &buffer, 1, MPI_INT, MPI_MAX, 0, splitcomm);
    PMPI_Comm_get_attr(splitcomm, namekey(), &communicator, &flag);

    MPI_Comm_rank(splitcomm,&rank);
    if ( flag ){
        /* relations between two calls, Bcast and Reduce, must be true */
        if ( rank !=0  ){
            assert(communicator->prim_bytes[Bcast] == communicator->prim_bytes[Reduce]);
            assert(communicator->prims[Bcast] == communicator->prims[Reduce]);
        }
    }

    MPI_Finalize();
    return 0;
}
