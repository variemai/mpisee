#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
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

    MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank / 2, &splitcomm);
    buffer = malloc(64);
    if ( rank  % 2 == 0  ){ //s1.0
        MPI_Send(buffer, 32, MPI_BYTE, rank+1, 0, MPI_COMM_WORLD); //Send to world 4x32
        MPI_Comm_rank(splitcomm, &rank);
        color = rank %2;
        if ( rank == 0 ){
            MPI_Sendrecv(buffer, 8, MPI_BYTE, rank+1, 0, buffer, 8, MPI_BYTE,
                         rank+1, 0, splitcomm, MPI_STATUS_IGNORE);
        }
        if ( rank == 1 ){
            MPI_Sendrecv(buffer, 8, MPI_BYTE, rank-1, 0, buffer, 8, MPI_BYTE,
                         rank-1, 0, splitcomm, MPI_STATUS_IGNORE);
        }
        MPI_Comm_split(splitcomm, color, rank / 2, &subcomm);
        MPI_Comm_rank(subcomm, &rank);
        MPI_Comm_free(&splitcomm);
        if ( rank < 1 && color == 0){
            MPI_Send(buffer,64,MPI_BYTE,rank+1,0,subcomm); //Send to s1.0_s2.0 1x64
        }
        else if ( rank > 1 && color == 1 ){
            MPI_Recv(buffer, 64, MPI_BYTE, rank-1, 0, subcomm, MPI_STATUS_IGNORE);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            printf("Process %d recvd %d bytes in subcomm color: %d\n",rank,64,color);
        }

    }
    else{
        MPI_Recv(buffer, 32, MPI_BYTE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Comm_rank(splitcomm, &rank);
        color = rank %2;
        MPI_Comm_split(splitcomm, color, rank / 2, &subcomm);
        MPI_Comm_rank(subcomm, &rank);
        if ( rank % 2 == 0 && color == 1){
            MPI_Send(buffer,16,MPI_BYTE,rank+1,0,subcomm); //Send to s1.1_s2.1 and 1x16
            MPI_Comm_free(&subcomm);
        }
        else if ( rank % 2 != 0 && color == 1 ){
            MPI_Recv(buffer, 16, MPI_BYTE, rank-1, 0,subcomm, MPI_STATUS_IGNORE);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            printf("Process %d recvd %d bytes in subcomm color: %d\n",rank,16,color);
        }
    }
    free(buffer);
    MPI_Finalize();
}
