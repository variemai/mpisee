#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#define NUMPROCS 8

int main (int argc, char *argv[]){
    int size;
    MPI_Comm splitcomm,subcomm;
    int rank,color;
    char *buffer;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if ( size != NUMPROCS ){
        if ( rank == 0 ){
            fprintf(stderr, "Error: Application runs only with %d processes\n",NUMPROCS);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    MPI_Comm_split(MPI_COMM_WORLD, rank % 4, rank / 4, &splitcomm);
    buffer = malloc(64);
    memset(buffer,0,sizeof(char)*64);
    MPI_Bcast(buffer, 64, MPI_CHAR, 0, splitcomm);
    MPI_Finalize();
    return 0;
}
