#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "../utils.h"
#include "../commprof.h"
#include <assert.h>

/* Metamorphic relations with Bcast function */

int main (int argc, char *argv[]){
    int flag,flag_N,rank;
    MPI_Comm splitcomm,splitcomm_N;
    int buffer[2];
    double buff;
    prof_attrs *communicator = NULL;
    prof_attrs *communicator_N = NULL;
    flag = 0;
    flag_N = 0;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    MPI_Comm_split(MPI_COMM_WORLD,rank % 2,rank,&splitcomm);
    buffer[0] = rank;
    buffer[1] = rank;
    /* 1st Bcast call with 2 MPI_INT */
    MPI_Bcast(&buffer, 2, MPI_INT, 0, splitcomm);
    PMPI_Comm_get_attr(splitcomm, namekey(), &communicator, &flag);


    MPI_Comm_split(MPI_COMM_WORLD,rank % 2,rank,&splitcomm_N);
    buff = (double) rank;
    /* 2nd Bcast call with 1 MPI DOUBLE and communicator */
    MPI_Bcast(&buff, 1, MPI_DOUBLE, 0, splitcomm_N);
    PMPI_Comm_get_attr(splitcomm, namekey(), &communicator_N, &flag_N);


    if ( flag && flag_N ){
        /* relation between the two communicators */
        assert(communicator->size == communicator_N->size);
        /* relations between two bcast calls must be true */
        assert(communicator->prim_bytes[Bcast] == communicator_N->prim_bytes[Bcast]);
        assert(communicator->prims[Bcast] == communicator_N->prims[Bcast]);
    }

    MPI_Finalize();
    return 0;
}
