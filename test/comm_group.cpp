#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#include <string.h>

int main (int argc, char *argv[]){
    int size;
    MPI_Comm newcomm,splitcomm,subcomm;
    MPI_Comm comm;
    MPI_Group group, newgroup, world_group;
    int rank, world_rank, world_size;
    char *buffer;
    const int ranks[2] = {0,1};
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank / 2, &splitcomm);
    buffer = malloc(64);
    if ( rank  % 2 == 0  ){
        MPI_Comm_group(splitcomm, &group);
        MPI_Comm_create(splitcomm, group, &subcomm);
        MPI_Comm_rank(subcomm, &rank);
        if ( rank % 2 == 0 ){
            MPI_Send(&size,1,MPI_INT,rank+1,0,subcomm);
        }
        else{
            MPI_Recv(&size, 1, MPI_INT, rank-1, 0, subcomm, MPI_STATUS_IGNORE);
        }
        if ( rank % 2 == 0 ){
            MPI_Send(buffer,64,MPI_BYTE,rank+1,0,subcomm);
        }
        else{
            MPI_Recv(buffer, 64, MPI_BYTE, rank-1, 0, subcomm, MPI_STATUS_IGNORE);
        }
    }
    else{
        MPI_Comm_group(splitcomm, &group);
        MPI_Comm_create(splitcomm, group, &subcomm);
        MPI_Comm_rank(subcomm, &rank);
        if ( rank % 2 == 0 ){
            MPI_Send(buffer,64,MPI_BYTE,rank+1,0,subcomm);
        }
        else{
            MPI_Recv(buffer, 64, MPI_BYTE, rank-1, 0, subcomm, MPI_STATUS_IGNORE);
        }
        MPI_Comm_split(subcomm, rank % 2, rank / 2, &newcomm);
        MPI_Comm_rank(newcomm, &rank);
        if ( rank % 2 == 0 ){
            MPI_Send(buffer,16,MPI_BYTE,rank+1,0,newcomm);
        }
        else{
            MPI_Recv(buffer, 16, MPI_BYTE, rank-1, 0, newcomm, MPI_STATUS_IGNORE);
        }

    }
    /* MPI_Comm_group(MPI_COMM_WORLD, &group); */
    /* MPI_Comm_size(MPI_COMM_WORLD, &world_size); */
    /* MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); */
    /* if ( MPI_Group_excl(world_group, 2, ranks, &newgroup) != MPI_SUCCESS ){ */
    /*     fprintf(stderr,"group_include_failed"); */
    /*     exit(EXIT_FAILURE); */
    /* } */
    /* if ( MPI_Group_incl(world_group, 2, ranks, &group) != MPI_SUCCESS ){ */
    /*     fprintf(stderr,"group_include_failed"); */
    /*     exit(EXIT_FAILURE); */
    /* } */
    /* if ( newgroup == MPI_GROUP_NULL ){ */
    /*         fprintf(stderr,"NULL NEWGROUP\n"); */
    /*         MPI_Finalize(); */
    /*         exit(EXIT_FAILURE); */
    /* } */
    /* if ( group == MPI_GROUP_NULL ){ */
    /*         fprintf(stderr,"NULL GROUP\n"); */
    /*         MPI_Finalize(); */
    /*         exit(EXIT_FAILURE); */
    /* } */
    /* // id = 0 */
    /* if ( world_rank != 0 && world_rank != 1 ){ */
    /*     if ( MPI_Comm_create_group(MPI_COMM_WORLD, newgroup, 0, &comm) != MPI_SUCCESS ){ */
    /*         fprintf(stderr,"comm_create_failed"); */
    /*         exit(EXIT_FAILURE); */
    /*     } */
    /*     MPI_Comm_size(comm,&size); */
    /*     MPI_Comm_rank(comm,&rank); */
    /*     printf("Communicator Y: World rank = %d, comm rank = %d, comm size = %d\n",world_rank,rank,size); */
    /*     if ( rank == 0 ){ */
    /*         MPI_Send(buffer,64,MPI_BYTE,rank+1,0,comm); */
    /*     } */
    /*     if ( rank == 1 ){ */
    /*         MPI_Recv(buffer, 64, MPI_BYTE, rank-1, 0, comm, MPI_STATUS_IGNORE); */
    /*     } */
    /* } */
    /* // id = 1 */
    /* if ( MPI_Comm_create(MPI_COMM_WORLD, group, &comm) != MPI_SUCCESS ){ */

    /*     fprintf(stderr,"comm_create_failed"); */
    /*     exit(EXIT_FAILURE); */
    /* } */
    /* if ( newcomm != MPI_COMM_NULL && newcomm != NULL){ */
    /*     MPI_Comm_rank(newcomm,&rank); */
    /*     MPI_Comm_size(newcomm,&size); */
    /*     if ( rank == 0 ){ */
    /*         MPI_Send(&size,1,MPI_INT,rank+1,0,newcomm); */
    /*     } */
    /*     if ( rank == 1 ){ */
    /*         MPI_Recv(&size, 1, MPI_INT, rank-1, 0, newcomm, MPI_STATUS_IGNORE); */
    /*     } */
    /*     printf("Communicator X: World rank = %d, newcomm rank = %d, newcomm size = %d\n",world_rank,rank,size); */
    /* } */
    /* if ( comm != MPI_COMM_NULL && comm != NULL){ */
    /*     MPI_Comm_rank(comm,&rank); */
    /*     MPI_Comm_size(comm,&size); */
    /*     printf("Communicator Y: World rank = %d, comm rank = %d, comm size = %d\n",world_rank,rank,size); */
    /*     if ( rank == 0 ){ */
    /*         MPI_Send(buffer,64,MPI_BYTE,rank+1,0,comm); */
    /*     } */
    /*     if ( rank == 1 ){ */
    /*         MPI_Recv(buffer, 64, MPI_BYTE, rank-1, 0, comm, MPI_STATUS_IGNORE); */
    /*     } */
    /* } */
    MPI_Finalize();
}

