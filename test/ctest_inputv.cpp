#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "../utils.h"

int random_int(int lower, int upper)
{
    int num = (rand() % (upper - lower + 1)) + lower;
    return num;
}


int main (int argc, char *argv[]){
    int i,j,len;
    int aac;
    char **aav = NULL;
    aac = random_int(1,MAX_ARGS);
    aav = (char**) malloc ( sizeof (char*)*aac);
    srand(time(NULL));
    for( i =0; i<aac; i++){
        aav[i] = (char*) malloc ( sizeof(char)*MAX_ARG_SIZE);
        len = random_int(0,MAX_ARG_SIZE-1);
        for ( j=0; j<len; j++ ){
            aav[i][j] = (char)random_int(33, 126);
        }
        aav[i][j]='\0';
    }
    MPI_Init(&aac, &aav);
    for ( i =0; i<aac; i++ ){
        assert ( strcmp(aav[i], av[i]) == 0 );
        printf("%s\n",av[i]);
    }

   MPI_Finalize();
   return 0;
}
