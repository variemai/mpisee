#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "../utils.h"
#include "mpi.h"

int main (int argc, char *argv[]){
    int i;
    MPI_Init(&argc, &argv);
    assert(ac > 0);
    for ( i =0; i<ac; i++ ){
        printf("AV[%d]: %s, ARGV[%d]: %s\n",i,av[i],i,argv[i]);
        //assert ( strcmp(av[i], argv[i]) == 0 );
    }

   MPI_Finalize();
   return 0;
}
