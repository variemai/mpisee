#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "../utils.h"
#include "../commprof.h"
#include <assert.h>

void check_names(prof_attrs *communicator, MPI_Comm comm){
    int bufflen,size,rank,i;
    char *buffer;
    char *recvbuffer;
    MPI_Comm_size(comm,&size);
    MPI_Comm_rank(comm,&rank);
    buffer = strdup(communicator->name);
    bufflen = strlen(buffer);
    recvbuffer = (char*) malloc (sizeof(char)*size*bufflen);
    MPI_Gather(buffer, bufflen, MPI_CHAR, recvbuffer, bufflen, MPI_CHAR, 0, comm);
    if ( rank == 0 ){
        for ( i = 0; i<size; i++ ){
            assert(strncmp(&recvbuffer[i*bufflen], buffer, bufflen) == 0 );
        }
    }
    free(buffer);
    free(recvbuffer);
}

/* Relations between MPI_Comm_create function calls */

int main (int argc, char *argv[]){
    int flag,rank;
    MPI_Comm comm0,comm1,comm2;
    MPI_Group commgrp;
    prof_attrs *communicator0 = NULL;
    prof_attrs *communicator1 = NULL;
    prof_attrs *communicator2 = NULL;
    flag = 0;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    MPI_Comm_group(MPI_COMM_WORLD,&commgrp);
    MPI_Comm_create(MPI_COMM_WORLD, commgrp, &comm0);
    MPI_Comm_create(MPI_COMM_WORLD, commgrp, &comm1);
    MPI_Comm_create(MPI_COMM_WORLD, commgrp, &comm2);

    PMPI_Comm_get_attr(comm0, namekey(), &communicator0, &flag);
    PMPI_Comm_get_attr(comm1, namekey(), &communicator1, &flag);
    PMPI_Comm_get_attr(comm2, namekey(), &communicator2, &flag);

    if ( flag ){
        assert(strcmp(communicator0->name,communicator1->name) != 0 );
        check_names(communicator0, MPI_COMM_WORLD);

        assert(strcmp(communicator0->name,communicator2->name) != 0 );
        check_names(communicator2, MPI_COMM_WORLD);

        assert(strcmp(communicator1->name,communicator2->name) != 0 );
        check_names(communicator1, MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return 0;
}
