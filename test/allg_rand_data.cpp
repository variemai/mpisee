#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "../utils.h"
#include "../commprof.h"
#include <assert.h>
#include <time.h>

int32_t maximum(int32_t a, int32_t b) { return((a) > (b) ? a : b); }

MPI_Datatype select_dt(int i){
    MPI_Datatype dt;
    switch (i) {
        case 0:
            dt = MPI_CHAR;
            break;
        case 1:
            dt = MPI_SHORT;
            break;
        case 2:
            dt = MPI_INT;
            break;
        case 3:
            dt = MPI_FLOAT;
            break;
        case 4:
            dt = MPI_DOUBLE;
            break;
        default:
            dt = MPI_CHAR;
            break;
    }
    return dt;
}
/* Metamorphic relations between Bcat and Reduce functions */

int main (int argc, char *argv[]){
    int rank,bytes,i,sz0,sz1,max,index[2];
    int cnt0,cnt1,size;
    int flag0,flag1;
    MPI_Comm comm0,comm1;
    char *sendbuf,*recvbuf;
    MPI_Datatype dt0,dt1;
    prof_attrs *communicator0 = NULL;
    prof_attrs *communicator1 = NULL;
    flag0 = 0;
    flag1 = 0;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    srand(time(NULL));
    if ( rank == 0 ){
        index[0] = rand() % 5;
        index[1] = rand() % 5;
        dt0 = select_dt(index[0]);
        dt1 = select_dt(index[1]);
        MPI_Type_size(dt0,&sz0);
        MPI_Type_size(dt1,&sz1);
        max = maximum(sz0,sz1);
        i = rand() % 15;
        bytes = 1;
        bytes = bytes<<i;
        while(bytes < max){
            i = rand() % 15;
            bytes = 1;
            bytes = bytes<<i;
        }

    }
    MPI_Bcast(index, 2, MPI_INT, 0, MPI_COMM_WORLD);
    if ( rank !=0  ){
        dt0 = select_dt(index[0]);
        dt1 = select_dt(index[1]);
        MPI_Type_size(dt0,&sz0);
        MPI_Type_size(dt1,&sz1);
    }
    MPI_Bcast(&bytes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    cnt0 = bytes/sz0;
    cnt1 = bytes/sz1;

    printf("cnt0=%d,sz0=%d,cnt1=%d,sz1=%d\n",cnt0,sz0,cnt1,sz1);
    sendbuf = (char*) malloc(bytes);
    MPI_Comm_split(MPI_COMM_WORLD,rank % 2,rank,&comm0);
    MPI_Comm_split(MPI_COMM_WORLD,rank % 2,rank,&comm1);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    recvbuf = (char*) malloc(bytes*size);
    memset(sendbuf,0,bytes);
    MPI_Allgather(sendbuf, cnt0, dt0, recvbuf, cnt0, dt0, comm0);
    MPI_Allgather(sendbuf, cnt1, dt1, recvbuf, cnt1, dt1, comm1);
    MPI_Scatter(sendbuf,cnt0,dt0,recvbuf,cnt0,dt0,0,comm0);
    MPI_Scatter(sendbuf,cnt1,dt1,recvbuf,cnt1,dt1,0,comm1);
    MPI_Gather(sendbuf,cnt0,dt0,recvbuf,cnt0,dt0,0,comm0);
    MPI_Gather(sendbuf,cnt1,dt1,recvbuf,cnt1,dt1,0,comm1);
    MPI_Bcast(sendbuf, cnt0, dt0, 0, comm0);
    MPI_Bcast(sendbuf, cnt1, dt1, 0, comm1);
    MPI_Allreduce(sendbuf, recvbuf, cnt0, dt0, MPI_MAX, comm0);
    MPI_Allreduce(sendbuf, recvbuf, cnt1, dt1, MPI_MAX, comm1);
    MPI_Reduce(sendbuf, recvbuf, cnt0, dt0, MPI_MAX, 0, comm0);
    MPI_Reduce(sendbuf, recvbuf, cnt1, dt1, MPI_MAX, 0, comm1);
    PMPI_Comm_get_attr(comm0, namekey(), &communicator0, &flag0);
    PMPI_Comm_get_attr(comm1, namekey(), &communicator1, &flag1);
    MPI_Comm_rank(comm0, &rank);
    if (flag1 && flag0) {

        assert(communicator0->prim_bytes[Allgather] == communicator1->prim_bytes[Allgather] &&
               communicator0->prim_bytes[Gather] == communicator1->prim_bytes[Gather] &&
               communicator0->prim_bytes[Scatter] == communicator1->prim_bytes[Scatter]);

        assert(communicator0->prim_bytes[Gather] == communicator0->prim_bytes[Allgather] &&
               communicator1->prim_bytes[Gather] == communicator1->prim_bytes[Allgather] &&
               communicator0->prim_bytes[Allreduce] == communicator0->prim_bytes[Reduce] &&
               communicator0->prim_bytes[Reduce] == communicator0->prim_bytes[Gather] &&
               communicator1->prim_bytes[Allreduce] == communicator1->prim_bytes[Reduce] &&
               communicator1->prim_bytes[Reduce] == communicator1->prim_bytes[Gather] );

        if ( rank != 0 ){
            assert(communicator0->prim_bytes[Gather] == communicator0->prim_bytes[Bcast] &&
                   communicator1->prim_bytes[Gather] == communicator1->prim_bytes[Bcast]);
        }
    }

    MPI_Finalize();
    return 0;
}
