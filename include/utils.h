/*****************************************************************************\
* This source file is part of mpisee
* Copyright (C) 2024  Ioannis Vardas - vardas@par.tuwien.ac.at
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <https://www.gnu.org/licenses/>
******************************************************************************/
#ifndef UTILS_H_
#define UTILS_H_

#include <cstdint>
#include <mpi.h>
#include <unistd.h>
#include <string>
#include <unordered_map>
#define NAMELEN 32
#define NUM_OF_PRIMS 76
#define MAX_ARG_STRING_SIZE 4096
#define MAX_ARGS 1024
#define MAX_ARG_SIZE 64
#define NUM_BUCKETS 8
#define MPISEE_MAJOR_VERSION 4
#define MPISEE_MINOR_VERSION 1
const int64_t buckets[NUM_BUCKETS-1] = {128,1024,8192,65536,262144,1048576,33554432};


// This is used to index the prim_names array
// Warning: This has to be in the same order as the prim_names array
enum primitives{
Send,
Recv,
Isend,
Irecv,
Sendrecv,
Isendrecv,
Ssend,
Issend,
Rsend,
Irsend,
Bsend,
Ibsend,
Waitall,
Wait,
Waitany,
Test,
Testany,
Testall,
Put,
Rput,
Get,
Rget,
Accumulate,
Raccumulate,
Fence,
Win_start,
Win_complete,
Win_post,
Win_wait,
Win_test,
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
Reduce_scatter_block,
Iallreduce,
Ibcast,
Ialltoall,
Iscatter,
Ibarrier,
Iallgather,
Iallgatherv,
Ialltoallv,
Ialltoallw,
Ireduce,
Igather,
Igatherv,
Iscan,
Iexscan,
Iscatterv,
Ireduce_scatter,
Ireduce_scatter_block,
Neighbor_allgather,
Neighbor_allgatherv,
Neighbor_alltoall,
Neighbor_alltoallv,
Neighbor_alltoallw,
Ineighbor_allgather,
Ineighbor_allgatherv,
Ineighbor_alltoall,
Ineighbor_alltoallv,
Ineighbor_alltoallw,
Init,
Init_thread
};

struct CommData {
    std::string name;
    int size;
};

struct DataEntry {
    int rank;
    int commId;
    int operationId;
    int bufferSizeMin;
    int bufferSizeMax;
    int calls;
    double time;
    uint64_t volume;
};

struct VolEntry {
    int operationId;
    int rank;
    int commId;
    uint64_t volume;
};

typedef struct profiler_attributes{
    char name[NAMELEN];
    int size;
    int comms;
    char id;
    double buckets_time[NUM_OF_PRIMS][NUM_BUCKETS];
    int buckets_msgs[NUM_OF_PRIMS][NUM_BUCKETS];
    uint64_t volume[NUM_OF_PRIMS][NUM_BUCKETS];
}prof_attrs;

typedef struct profiler_data{
    char name[NAMELEN];
    int size;
    double buckets_time[NUM_OF_PRIMS][NUM_BUCKETS];
    int buckets_msgs[NUM_OF_PRIMS][NUM_BUCKETS];
    uint64_t volume[NUM_OF_PRIMS][NUM_BUCKETS];
}prof_data;

typedef struct request_data{
    MPI_Request *req;
    MPI_Comm comm;
}rq;


// Metadata associated with an MPI communicator
typedef struct profiler_metadata {
    char name[NAMELEN];
    int size;
    int comms;
    char id;
} prof_metadata;

// Structure to store profiling data for a primitive and bucket
typedef struct PrimBucketInfo {
    double time;
    int num_messages;
    uint64_t volume;
}primBucketInfo;

typedef struct comm_profiler {
    std::unordered_map<int, primBucketInfo> data_map;
} comm_profiler;


typedef struct comm_data{
    int comm_id;
    int prim;
    int bucketIndex;
    double time;
    int num_messages;
    uint64_t volume;
} comm_data;

typedef struct comm_all{
    char name[NAMELEN];
    int size;
    int datasize;
} comm_all;

typedef struct comm_profiler_meta_pair{
    comm_profiler prof;
    prof_metadata meta;
} prof_meta_pair;


/* Globals */
extern char *appname;
extern const char prim_names[][NUM_OF_PRIMS];
extern int ac;
extern char *av[MAX_ARGS];

void mcpt_abort (const char *fmt, ...);

void getProcCmdLine (int *ac, char **av);
char * get_appname(void);

void getRunCmd(int argc, char **argv);

#endif // UTILS_H_
