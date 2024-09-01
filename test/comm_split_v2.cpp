#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#define NUMPROCS 16

int main (int argc, char *argv[]){
    int size;
    MPI_Comm splitcomm,subcomm,leafcomm;
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

    MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank / 2, &splitcomm); //s0
    buffer = malloc(64);
    if ( rank  % 2 == 0  ){ //s1_0
        MPI_Send(buffer, 32, MPI_BYTE, rank+1, 0, MPI_COMM_WORLD);
        MPI_Comm_rank(splitcomm, &rank);
        color = rank %2;
        if ( rank ==  0 ){
            MPI_Send(buffer,64,MPI_BYTE,rank+1,0,splitcomm);
        }
        if ( rank == 1 ){
            MPI_Recv(buffer, 64, MPI_BYTE, rank-1, 0, splitcomm, MPI_STATUS_IGNORE);
        }

        MPI_Comm_split(splitcomm, color, rank / 2, &subcomm);//s1.0_s2.0
        MPI_Comm_rank(subcomm, &rank);
        if ( color == 0 ){ //s1_0 parent: s0_0
            /* MPI_Comm_free(&splitcomm); */
            if ( rank == 0 ){
                MPI_Send(buffer,8,MPI_BYTE,rank+1,0,subcomm);
            }
            if(rank == 1){
                MPI_Recv(buffer, 8, MPI_BYTE, rank-1, 0, subcomm, MPI_STATUS_IGNORE);
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                printf("Process %d recvd %d bytes in subcomm color: %d\n",rank,8,color);
            }
        }
        else{ //s1.0_s2.1
            color = rank % 2;
            MPI_Comm_split(subcomm, color, rank / 2, &leafcomm);
            MPI_Comm_rank(leafcomm, &rank);
            if ( rank < 1 ) //s1.0_s2.1_s3.0 and s1.0_s2.1_s3.1
                MPI_Send(buffer,4,MPI_BYTE,rank+1,0,leafcomm);
            else
                MPI_Recv(buffer, 4, MPI_BYTE, rank-1, 0, leafcomm, MPI_STATUS_IGNORE);
        }

    }
    else{ //s1_1
        MPI_Recv(buffer, 32, MPI_BYTE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Comm_rank(splitcomm, &rank);
        if  ( rank % 2 == 0  ){
            MPI_Send(buffer,16,MPI_BYTE,rank+1,0,splitcomm);
        }
        else{
            MPI_Recv(buffer, 16, MPI_BYTE, rank-1, 0,splitcomm, MPI_STATUS_IGNORE);
        }
        color = rank %2;
        MPI_Comm_split(splitcomm, color, rank / 2, &subcomm);
        MPI_Comm_rank(subcomm, &rank);
        if ( rank % 2 == 0 ){
            MPI_Send(buffer,4,MPI_BYTE,rank+1,0,subcomm);
            /* MPI_Comm_free(&subcomm); */
        }
        else{
            MPI_Recv(buffer, 4, MPI_BYTE, rank-1, 0,subcomm, MPI_STATUS_IGNORE);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            printf("Process %d recvd %d bytes in subcomm color: %d\n",rank,4,color);
        }
        if ( color == 1 ){
            color = rank % 2;
            MPI_Comm_split(subcomm, color, rank / 2, &leafcomm); //s2
            MPI_Comm_rank(leafcomm, &rank);
            if ( color == 1 ){
                if ( rank < 1 )
                    MPI_Send(buffer,8,MPI_BYTE,rank+1,0,leafcomm);
                else
                    MPI_Recv(buffer, 8, MPI_BYTE, rank-1, 0, leafcomm, MPI_STATUS_IGNORE);
            }
        }
    }
    free(buffer);
    MPI_Finalize();
}

