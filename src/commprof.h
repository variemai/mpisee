#ifndef COMMPROF_H_
#define COMMPROF_H_

#include "mpi.h"
#include <unordered_map>
#define MAX_DIMS 8
#include "utils.h"
#include <vector>
extern int prof_enabled;
extern prof_attrs **local_data;
// extern prof_attrs **local_comms;
extern std::vector<MPI_Comm> comms_table;
extern std::unordered_map<MPI_Request, MPI_Comm> requests_map;
extern std::unordered_map<MPI_Win, MPI_Comm> comm_map;
extern std::vector<std::pair<prof_meta_pair*, MPI_Group>> free_array;

extern int local_cid;
extern int my_coms;
extern int keys[2]; // keyval[0]  contains metadata
                      // keyval[1]  contains profiling data

extern "C" {
int
namekey(void);
}
extern "C" {
int
namedel(MPI_Comm comm, int keyval, void *attr, void *s);
}

extern "C" {
void
profile_this(MPI_Comm comm, int64_t count,MPI_Datatype datatype,int prim,
             double t_elapsed,int v);
}

extern "C" {
int
win_namekey(void);
}
#endif
