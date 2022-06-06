#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "../utils.h"

int main (int argc, char *argv[]){
    int i;
    MPI_Init(&argc, &argv);
    getRunCmd(argc, argv);
    assert(ac > 0);
    for ( i =0; i<ac; i++ ){
        assert ( strcmp(av[i], argv[i]) == 0 );
        printf("%s\n",av[i]);
    }

  for ( i =0; i<ac; i++ ){
      free(av[i]);
  }
  free(av);
   MPI_Finalize();
    return 0;
}
