#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#include <string.h>
#define NUMPROCS 16

static int namedel(MPI_Comm comm, int keyval, void *attr, void *s) {
  char *com = (char*)attr;
  if ( comm == NULL )
      return -1;
  else
      free(com);
  return 0;
}

static int namekey() {
  // hidden key value for type attributes
  static int namekeyval = MPI_KEYVAL_INVALID;

  if (namekeyval == MPI_KEYVAL_INVALID) {
    MPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN,namedel,&namekeyval,NULL);
  }

  return namekeyval;
}

int main (int argc, char *argv[]){
    int size;
    MPI_Comm subcomm,dupcomm;
    MPI_Group group;
    int rank,flag;
    char *buffer,*buf,*b;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if ( size != NUMPROCS ){
        if ( rank == 0 ){
            fprintf(stderr, "Error: Application runs only with %d processes\n",NUMPROCS);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    MPI_Comm_group(MPI_COMM_WORLD, &group);
    MPI_Comm_create(MPI_COMM_WORLD, group, &subcomm);
    buffer = malloc(64);
    strcpy(buffer, "PAPARIA");
    MPI_Comm_set_attr(subcomm, namekey(), buffer);
    MPI_Comm_dup(subcomm, &dupcomm);
    b = malloc(64);
    strcpy(b, "TROLL");
    MPI_Comm_set_attr(dupcomm, namekey(), b);
    if ( rank == 0 ){
        MPI_Comm_get_attr(subcomm, namekey(), &buf, &flag);
        if ( flag )
            printf("Subcomm: %s\n",buf);
        MPI_Comm_get_attr(dupcomm, namekey(), &buf, &flag);
        if ( flag )
            printf("Dupcomm: %s\n",buf);
    }
    free(buffer);
    MPI_Finalize();
}
