#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "../utils.h"
#include "../commprof.h"
#include <assert.h>

int ctest_ac;
char **ctest_av;

int main (int argc, char *argv[]){
    int size,i,rank;
    char *ptr = NULL;
    MPI_Init(&argc, &argv);
    ctest_av = (char**) malloc ( sizeof(char)*argc*64 );
    ptr = ctest_av[0];
    ctest_ac = argc;
    for ( i =0; i<argc; i++ ){
      ctest_av[i]=strdup(argv[i]);
    }
    for ( i =0; i<argc; i++ ){
      printf("%s ",ctest_av[i]);
    }
    MPI_Finalize();
    return 0;
}
