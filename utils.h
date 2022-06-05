#ifndef UTILS_H_
#define UTILS_H_

#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#define NAMELEN 256
#define PRIMLEN 16
#define NUM_OF_PRIMS 36

enum primitives{
Send,           /* DO NOT ADD ANYHTING BEFORE THIS */
Recv,
Isend,
Irecv,
Ssend,
Issend,
Probe,
Iprobe,
Waitall,
Wait,
Waitany,
Test,
Testany,
Testall,
Sendrecv,
Bcast,
Barrier,
Allreduce,
Allgather,
Allgatherv,
Alltoall,
Alltoallv,
Alltoallw,
Reduce,
Gather,
Gatherv,
Scan,
Exscan,
Scatter,
Scatterv,
Reduce_scatter,
Iallreduce,
Ibcast,
Ialltoall,
Iscatter,
Ibarrier   /* DO NOT ADD ANYTHING AFTER THIS */
};


typedef struct profiler_attributes{
    char name[NAMELEN];
    uint64_t bytes;
    uint64_t msgs;
    int size;
    uint32_t prims[NUM_OF_PRIMS];
    uint64_t prim_bytes[NUM_OF_PRIMS];
    double time_info[NUM_OF_PRIMS];
}prof_attrs;

typedef struct request_data{
    MPI_Request *req;
    MPI_Comm comm;
}rq;

/* Globals */
extern char *appname;
extern const char prim_names[][NUM_OF_PRIMS];

void mcpt_abort (char *fmt, ...);

void getProcCmdLine (int *ac, char **av);
char * get_appname(void);

#endif // UTILS_H_
