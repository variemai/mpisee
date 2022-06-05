#include <stdio.h>
#include <assert.h>
#include "arith.h"
#include <stdlib.h>
#include <string.h>
#include "except.h"
#include "mem.h"
#include "stack.h"
#include "queue.h"
#include "list.h"
#include "atom.h"
#include <time.h>
#include "table.h"
#include <mpi.h>
#include "../utils.h"


Except_T FileOpen_failed = {"FILE OPEN FAILED"};

typedef struct _node{
    char *name;
    struct _node *next;
}node;



typedef struct fuck_my_life{
    const char* key;
    int* value;
}fml;

node* head;

void newNode(const char *n){
    node *ptr;
    ptr = ALLOC(sizeof(node));
    ptr->name = ALLOC(80);
    ptr->next = head->next;
    head->next = ptr;
    ptr->name=strcpy(ptr->name,n);
}

void applyPrint(void **ptr,void *cl){
    printf("%s\n",(char*)*ptr);
}

void apply_print(const void *key,void **value,void *cl){
    int *val = *value;
    printf("%s %d\n",(char*)key,*val);
}

int compare(const void *x, const void *y) {
    return strcmp(*(char **)x, *(char **)y);
}

int main(int32_t argc, char const *argv[]){
	/* uint32_t i; */
	/* const char *strn; */
	/* clock_t begin,end; */
	/* double time_spent; */
    /* char **buf; */
    int length,i;
    const char *com0, *com1, *com2, *com3;
    /* int* b_com0, *b_com1, *b_com2; */
    int *b_com3;
    void** array;
    /* List_T lista; */
    Table_T table;
    MPI_Comm newcomm,split,third;
    MPI_Request req0, req1, req2;
    /* prof_attrs *com_0, *com_1, *com_2; */
    int rank;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    table = Table_new(5,NULL,NULL);
    MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank, &newcomm);
    /* com_0 = (prof_attrs*) malloc ( sizeof(prof_attrs) ); */
    MPI_Comm_rank(newcomm, &rank);
    /* NEW(b_com0); */
    /* NEW(b_com1); */
    /* NEW(b_com2); */
    com0=Atom_string("com0");
    com1=Atom_string("com1");
    com2=Atom_string("com2");
    com3=Atom_string("com3");
    NEW(b_com3);
    if ( rank == 0 ){
        MPI_Isend(&rank, 1, MPI_INT, 1, 0, newcomm, &req0);
        Table_put(table, &req0, (char*)com0);
    }
    if ( rank == 1 ){
        MPI_Irecv(b_com3, 1, MPI_INT, 0, 0, newcomm, &req0);
        Table_put(table, &req0, (char*)com0);
    }
    /* array = Table_toArray(table, NULL); */
    /* for (i = 0; array[i]; i += 2){ */
    /*     printf("%d: %p\t%p\n", rank,(MPI_Request)array[i],(MPI_Comm)array[i+1]); */
    /*     fflush(stdout); */
    /* } */
    /* FREE(array); */
    MPI_Wait(&req0, MPI_STATUS_IGNORE);
    /* *b_com0 = 100; */
    /* *b_com1 = 10; */
    /* *b_com2 = 33; */
    *b_com3 = 20;
    MPI_Comm_split(MPI_COMM_WORLD, rank %2, rank, &split);
    MPI_Comm_rank(split, &rank);
    if ( rank == 0 ){
        MPI_Isend(&rank, 1, MPI_INT, 1, 0, split, &req1);
        Table_put(table, &req1, (char*)com2);
    }
    if ( rank == 1 ){
        MPI_Irecv(b_com3, 1, MPI_INT, 0, 0, split, &req1);
        Table_put(table, &req1, (char*)com2);
    }
    MPI_Wait(&req1, MPI_STATUS_IGNORE);
    req2 = MPI_REQUEST_NULL;
    MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank, &third);
    Table_put(table, &req2, (char*)com3);
    Table_remove(table, &req1);
    length = Table_length(table);
    printf("Table length = %d\n",Table_length(table));
    array = Table_toArray(table, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for (i = 0; array[i]; i += 2){
        printf("%d: %p\t%s\n", rank,(MPI_Request)array[i],(char*)array[i+1]);
        fflush(stdout);
    }
    /* qsort(array, Table_length(table), 2*sizeof (*array), */
    /*    compare); */
    /* if ( rank == 0 ) */
    /*   printf("%\n", MPI_REQUEST_NULL); */
    MPI_Barrier(MPI_COMM_WORLD);
   FREE(array);
   /* FREE(b_com0); */
   /* FREE(b_com1); */
   /* FREE(b_com2); */
   FREE(b_com3);
   MPI_Finalize();
   return 0;

}
