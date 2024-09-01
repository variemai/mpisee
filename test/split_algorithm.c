#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <ctype.h>
#include <unistd.h>
#define RANKS 6

static int
compare_int(const void *a, const void *b){
    int *a0, *b0;
    a0 = (int*)a;
    b0 = (int*)b;
    return (*a0)-(*b0);
}

int main(int argc,char *argv[]){
    int color[6] = {19,30,999,19,999,40};
    int size,rank,*recvbuf, *uniq,tmp,j,comms,i;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    recvbuf = (int*) malloc ( sizeof(int)*size );
    uniq = (int*) malloc ( sizeof(int)*size );
    MPI_Allgather(&color[rank], 1, MPI_INT, recvbuf, 1, MPI_INT, MPI_COMM_WORLD);
    qsort(recvbuf, size, sizeof(int), compare_int);
    j = 0;
    tmp = 0;
    comms = 0;
    for ( i = 0; i<size; i++ ){
        if ( recvbuf[i] > tmp  ){
            comms ++ ;
            tmp = recvbuf[i];
            uniq[j] = tmp;
            j++;
        }
    }
    for ( j =0; j<comms; j++ ){
        if ( color[rank] == uniq[j]  )
            break;
    }
    printf ( "Rank %d my color = %d\n",rank,uniq[j] );
    if ( rank == 0 ){
        for ( i = 0; i < comms; i++ ){
            printf("Color %d\n",uniq[i]);
        }
        printf("Total unique colors %d\n",comms);
    }
    free(uniq);
    free(recvbuf);
    MPI_Finalize();

    return 0;
}
