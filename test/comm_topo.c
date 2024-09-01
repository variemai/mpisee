/* virt_topo.c
* https://computing.llnl.gov/tutorials/mpi/#Virtual_Topologies
*/
#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#include <string.h>
#define SIZE 16
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

int main (int argc, char *argv[]){
    int rank = 0, source, dest, outbuf, i;
    MPI_Comm cartcomm,splitcomm,dupcomm0,dupcomm1;
    int world_rank, world_size;
    int periods[2] = {1,1};
    int dims[2] = {4,2};
    MPI_Request reqs[4];
    MPI_Status stats[4];
    int inbuf[4]={MPI_PROC_NULL,MPI_PROC_NULL, MPI_PROC_NULL,MPI_PROC_NULL},
        nbrs[4], reorder=0, coords[2];

    int tag = 1;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_split(MPI_COMM_WORLD, world_rank % 2, world_rank / 2, &splitcomm);
    memset(nbrs, 0, sizeof(int)*4);
    if ( world_size == SIZE ){
        MPI_Cart_create(splitcomm, 2, dims, periods,reorder, &cartcomm);
        MPI_Comm_rank(cartcomm, &rank);
        MPI_Cart_coords(cartcomm, rank, 2, coords);
        MPI_Cart_shift(cartcomm, 0, 1, &nbrs[UP], &nbrs[DOWN]);
        /* MPI_Cart_shift(cartcomm, 1, 1, &nbrs[LEFT], &nbrs[RIGHT]); */
        /* printf("P[%d]: coords= %d %d \n", */
        /*        rank,coords[0],coords[1]); */
        /* printf("P[%d]: neighbors(u,d,l,r)=%d %d %d %d\n", */
               /* rank,nbrs[UP],nbrs[DOWN],nbrs[LEFT], nbrs[RIGHT]); */
        if ( world_rank % 2 == 0 ){
            MPI_Comm_dup(splitcomm, &dupcomm0);
            outbuf = rank;
            for (i=0; i<2; i++) {
                dest = nbrs[i];
                source = nbrs[i];
                MPI_Isend(&outbuf, 1, MPI_INT, dest, tag,
                          cartcomm, &reqs[i]);
                MPI_Irecv(&inbuf[i], 1, MPI_INT, source, tag,
                          cartcomm, &reqs[i+2]);
            }
            MPI_Waitall(4, reqs, stats);
        }
        else{
            outbuf = rank*2;
            for (i=0; i<2; i++) {
                dest = nbrs[i];
                source = nbrs[i];
                MPI_Isend(&outbuf, 1, MPI_INT, dest, tag,
                          cartcomm, &reqs[i]);
                MPI_Irecv(&inbuf[i], 1, MPI_INT, source, tag,
                          cartcomm, &reqs[i+2]);
            }
            MPI_Waitall(4, reqs, stats);
        }
        printf("rank= %d inbuf(u,d,l,r)= %d %d %d %d\n",
               world_rank,inbuf[UP],inbuf[DOWN],inbuf[LEFT],inbuf[RIGHT]);
    }
    else{
        printf("Must specify %d processors. Terminating.\n",SIZE);
    }
    MPI_Comm_dup(MPI_COMM_WORLD, &dupcomm1);

    MPI_Finalize();
}

