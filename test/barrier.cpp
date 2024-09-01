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
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
