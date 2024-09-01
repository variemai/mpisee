#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#include <string.h>

typedef struct comm_p_info{
    char *name;
    unsigned long long bytes;
    char *creator;
}cpi;

static int namedel(MPI_Comm comm, int keyval, void *attr, void *s) {
    cpi *com = (cpi*)attr;
    free(com->name);
    free(com->creator);
    free(com);
    return MPI_SUCCESS;
}

static int namekey() {
  // hidden key value for type attributes
  static int namekeyval = MPI_KEYVAL_INVALID;

  if (namekeyval == MPI_KEYVAL_INVALID) {
    MPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN,namedel,&namekeyval,NULL);
  }

  return namekeyval;
}



int comm_create_wrapper(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm){
    int ret;
    char *buf, *buffer = NULL;
    int len;
    cpi *com_info;
    ret = PMPI_Comm_create(comm, group, newcomm);
    //check here
    com_info = malloc (sizeof(cpi));
    com_info->bytes = 0;
    com_info->name = malloc(32);
    com_info->creator = malloc(32);
    strcpy(com_info->name, "COMM_0");
    strcpy(com_info->creator, "COMM_CREATE");
    PMPI_Comm_set_attr(*newcomm,namekey(),com_info);
    return ret;
}

int send_wrapper(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm){
    /* MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) */
    cpi *com_info;
    int ret,flag,size;
    ret = PMPI_Send(buf, count,  datatype,  dest,  tag,  comm);
    PMPI_Comm_get_attr(comm, namekey(), &com_info, &flag);
    if ( flag ){
        printf("%s HAD %llu\n",(*com_info).name,(*com_info).bytes);
        MPI_Type_size(datatype, &size);
        com_info->bytes += (count*size);
    }
    return ret;
}


void another(MPI_Comm newcome){
    printf("%p\n",newcome);
}

int main (int argc, char *argv[]){
    int size;
    char *comm_name;
    cpi *com_info;
    int len;
    MPI_Comm newcomm,splitcomm;
    MPI_Group group;
    int flag,key;
    int rank,world_rank;
    char *buffer;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    /* printf("Hello from all %d ranks\n",size); */
    /* MPI_Comm_get_name(MPI_COMM_WORLD, comm_name, &len); */
    MPI_Comm_group(MPI_COMM_WORLD, &group);
    /* wrapper_shit(MPI_COMM_WORLD, group, &newcomm); */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    /* MPI_Comm_create_group(MPI_Comm comm, MPI_Group group, int tag, MPI_Comm *newcomm) */
    /* if ( rank == 0 || rank == 1 ){ */
    MPI_Comm_create(MPI_COMM_WORLD, group,&newcomm);
    MPI_Comm_rank(newcomm, &rank);
    /* } */
    /* MPI_Comm_set_name(newcomm, "BOURDELO"); */
    /* MPI_Comm_get_name(newcomm, comm_name, &len); */
    /* set_attr_wrapper(newcomm, 0, val); */
    /* MPI_Comm_set_attr(newcomm, namekey(), com_info); */
    /* comm_create_wrapper(MPI_COMM_WORLD, group, &newcomm); */
    /* bullshit(&newcomm) ; */
    /* sprintf(comm_name, "TRIBOURDELO"); */
    /* MPI_Comm_set_attr(newcomm, 0, comm_name); */
    MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank / 2, &splitcomm);
    buffer = malloc(64);
    if ( rank == 0  ){
        /* MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); */
        /* printf("Rank = %d, World = %d\n",rank,world_rank); */
        MPI_Send(&size,1,MPI_INT,1,0,newcomm);
        MPI_Recv(buffer, 64, MPI_CHAR, 1, 1, newcomm, MPI_STATUS_IGNORE);
        MPI_Send(&size,1,MPI_INT,1,0,splitcomm);
        /* printf("I received %s from rank %d\n",buffer,rank); */
        /* MPI_Send(&size,1,MPI_INT,1,2,MPI_COMM_WORLD); */
        /* send_wrapper(&size, 1, MPI_INT, 1, 0, newcomm); */
        /* MPI_Comm_get_attr(newcomm, namekey(), &com_info, &flag); */
        /* if ( flag ){ */
        /*     printf("Rank %d: %s HAS %llu\n",rank,com_info->name,com_info->bytes); */
        /* } */
    }
    if ( rank == 1 ){
        MPI_Recv(&size, 1, MPI_INT, 0, 0, newcomm, MPI_STATUS_IGNORE);
        /* strcpy(buffer, "This is a shitty application"); */
        MPI_Send(buffer,64,MPI_CHAR,0,1,newcomm);
        /* MPI_Recv(&size, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE); */
    }
    if ( rank == 2 ){
        /* strcpy(buffer, "This is a shitty application"); */
        MPI_Send(buffer,64,MPI_CHAR,3,3,MPI_COMM_WORLD);
        MPI_Recv(&size, 1, MPI_INT, 0, 0, splitcomm, MPI_STATUS_IGNORE);
        printf("I received size = %d from rank 0\n",size);
    }
    if ( rank == 3 ){
        MPI_Recv(buffer, 64, MPI_CHAR, 2, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        /* printf("I received %s from rank 2\n",buffer); */
    }
    /* if ( rank ==0  ){ */
    /*     MPI_Comm_get_name(*newcomm, buffer, &len); */
    /*     printf("Name = %.*s\n",len,comm_name); */
    /* } */
    /*     printf("WORLD = %ld, BOURDELO = %ld\n",MPI_COMM_WORLD,newcomm); */
    /* } */
    MPI_Finalize();
}
