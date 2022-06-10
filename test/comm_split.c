#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "../utils.h"
#include "../commprof.h"
#include <assert.h>


/* Relations between MPI_Comm_split function calls */

int main (int argc, char *argv[]){
    int flag,rank;
    MPI_Comm comm0,comm1,comm2;
    prof_attrs *communicator0 = NULL;
    prof_attrs *communicator1 = NULL;
    prof_attrs *communicator2 = NULL;
    flag = 0;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank, &comm0);
    MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank, &comm1);
    MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank, &comm2);

    PMPI_Comm_get_attr(comm0, namekey(), &communicator0, &flag);
    PMPI_Comm_get_attr(comm1, namekey(), &communicator1, &flag);
    PMPI_Comm_get_attr(comm2, namekey(), &communicator2, &flag);

    if ( flag ){
        assert(strcmp(communicator0->name,communicator1->name) != 0 );

        assert(strcmp(communicator0->name,communicator2->name) != 0 );

        assert(strcmp(communicator1->name,communicator2->name) != 0 );
    }
    MPI_Finalize();
    return 0;
}
