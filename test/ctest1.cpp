#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "../utils.h"
#include "../commprof.h"
#include <assert.h>

int main (int argc, char *argv[]){
    int flag,rank,recv;
    MPI_Comm splitcomm;
    prof_attrs *communicator = NULL;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_split(MPI_COMM_WORLD,rank % 2,rank,&splitcomm);
    assert(my_coms == 2);
    MPI_Allreduce(&rank,&recv,1,MPI_INT,MPI_MIN,splitcomm);
    flag = 0;
    PMPI_Comm_get_attr(splitcomm, namekey(), &communicator, &flag);
    if ( flag ){
        assert(communicator->prims[Allreduce] == 1);
        assert(communicator->time_info[Allreduce] > 0);
        assert(communicator->prim_bytes[Allreduce] == sizeof(int));
    }
    MPI_Finalize();
    return 0;
}
