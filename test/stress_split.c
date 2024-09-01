#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#define ITERS 512              /* Default number of iterations */

int main (int argc, char *argv[]){
    int size, rank,iters,i;
    MPI_Comm splitcomm;
    MPI_Init(&argc,&argv);
    if ( argc < 2 ){
        iters = ITERS;
    }
    else{
        iters = atoi(argv[1]);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for ( i = 0; i<iters; i++ ){
        MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank / 2, &splitcomm);
        MPI_Comm_free(&splitcomm);
    }
    MPI_Finalize();
    return 0;
}
