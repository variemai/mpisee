#include "utils.h"
#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <mpi.h>
#include <ostream>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <cstdarg>
#include <unordered_map>
#include <vector>
#include "commprof.h"
#include <algorithm>
#include "symbols.h"
#include "create_db.h"
#include "utils.h.in"
#include <iostream>
#include <chrono>
#include "mpisee_fortran.h"

// int global_rank; // For debugging purposes
int prof_enabled = 1;
int local_cid= 0;
int my_coms = 0;
int ac;
char *av[MAX_ARGS];
int keys[2]; // keyval[0]  contains metadata
               // keyval[1]  contains profiling data
               
/* Necessary bookkeeping data structures */
std::unordered_map<MPI_Request, MPI_Comm> requests_map;
std::unordered_map<MPI_Win, MPI_Comm> comm_map;
std::vector<MPI_Comm> comms_table;
std::vector<std::pair<MPI_Group, prof_attrs*>> group_table;
std::vector<std::pair<prof_meta_pair*, MPI_Group>> free_array;


/* Tool date */
char mpisee_build_date[sizeof(__DATE__)] = __DATE__;
char mpisee_build_time[sizeof(__TIME__)] = __TIME__;
double total_time = 0.0;
double init_time;

extern "C" {
int
namedel(MPI_Comm comm, int keyval, void *attr, void *s)
{
  prof_attrs *com = (prof_attrs*)attr;
  free(com);
  return MPI_SUCCESS;
}
}

extern "C" {
int
namekey(void)
{
  // hidden key value for type attributes
  static int namekeyval = MPI_KEYVAL_INVALID;

  if (namekeyval == MPI_KEYVAL_INVALID) {
    PMPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN,namedel,&namekeyval,NULL);
  }

  return namekeyval;
}
}

extern "C" {
int
win_namedel(MPI_Comm comm, int keyval, void *attr, void *s){
    MPI_Comm *com = (MPI_Comm*)attr;
    free(com);
    return MPI_SUCCESS;
}
}


extern "C" {
int
win_namekey(void){
    static int win_keyval = MPI_KEYVAL_INVALID;

    if (win_keyval == MPI_KEYVAL_INVALID) {
#ifdef MPICH_NAME
        MPI_Win_create_keyval(MPI_WIN_NULL_COPY_FN, win_namedel, &win_keyval, NULL);
#endif
#ifdef OMPI_MAJOR_VERSION
        PMPI_Win_create_keyval(MPI_WIN_NULL_COPY_FN, MPI_WIN_NULL_DELETE_FN, &win_keyval, NULL);
#endif

    }
    return win_keyval;

}
}

int getPrimBucketKey(int prim, int bucketIndex)
{
    //std::cout << "mpisee: prim = " << prim << ",
    // bucketIndex = " << bucketIndex << " key = " << prim * NUM_BUCKETS + bucketIndex << std::endl;
    return prim * NUM_BUCKETS + bucketIndex;
}

// Initialize the profiling structure for the communicator
// Called by communicator creation functions
void
alloc_init_commprof(MPI_Comm comm, char c)
{

    prof_metadata *metadata;
    comm_profiler *prof;
    metadata = new prof_metadata();
    prof = new comm_profiler();
    PMPI_Comm_size(comm, &metadata->size);
    metadata->comms = local_cid++;
    metadata->id = c;
    PMPI_Comm_set_attr(comm, keys[0], metadata);
    PMPI_Comm_set_attr(comm, keys[1], prof);

}


int
choose_bucket(int64_t bytes) {
    int index;
    for (index = 0; index < NUM_BUCKETS-1; index++) {
        if ( buckets[index] > bytes) {
            return index;
        }
    }
    return NUM_BUCKETS-1;
}


void insertOrUpdatePrimBucketInfo(std::unordered_map<int, primBucketInfo>& map,
                                  int key, double time, uint64_t volume) {

    // Create a new primBucketInfo object
    primBucketInfo newInfo;
    newInfo.time = time;
    newInfo.num_messages = 1;
    newInfo.volume = volume;

    // Check if the key exists in the map
    auto it = map.find(key);

    if (map.count(key) > 1) {
        std::cerr << "Hash collision detected for key " << key << std::endl;
    }

    if (it == map.end()) {
        // Key not found, insert the new pair
        map[key] = newInfo;
        //std::cout << "mpisee: Inserted new key " << key << std::endl;
    } else {
        // Key found, update the existing value
        it->second.time += newInfo.time;
        it->second.num_messages += newInfo.num_messages;
        it->second.volume += newInfo.volume;
        //std::cout << "mpisee: Updated key " << key << std::endl;
    }
}

// Profile the communication
extern "C" {
void
profile_this(MPI_Comm comm, int64_t count,MPI_Datatype datatype,int prim,
             double t_elapsed,int v){


    int size,flag,bucket_index;
    int64_t sum;
    comm_profiler *comm_prof;
    PMPI_Comm_get_attr(comm, keys[1], &comm_prof, &flag);
    /* Debugging code
    int rank;
    PMPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (flag) {
        printf("Rank %d: Found the map\n", rank);
    } else {
        printf("Rank %d: Map not found\n", rank);
    }
    printf("Rank %d: Inserting data into map\n", rank);

    End of debugging code */

    if ( datatype != MPI_DATATYPE_NULL ){
        PMPI_Type_size(datatype, &size);
        sum = count * size;
    }
    else{
        sum = count;
    }
    if (flag) {
        if ( v == 0 ){
            bucket_index = choose_bucket(sum);
            insertOrUpdatePrimBucketInfo(comm_prof->map,
                                         getPrimBucketKey(prim, bucket_index),
                                         t_elapsed,  sum);
        }
        // Don't record the buffer range for [v,w] collectives
        else{
            insertOrUpdatePrimBucketInfo(comm_prof->map,
                                         getPrimBucketKey(prim, 0),
                                         t_elapsed, sum);
        }
    }
    else{
        mcpt_abort("empty flag when profiling %s - this might be a bug\n",prim_names[prim]);
    }
}
}

// Compose the communicator's name
// Called only in Finalize and Comm_free
void
overwrite_name(prof_metadata **com_prof, int id0, int id1) {
    int requiredSize,r;
    requiredSize = snprintf(NULL,0, "%c%d.%d",(*com_prof)->id, id0, id1);
    if (requiredSize < 0 || requiredSize + 1 >= NAMELEN) {
        mcpt_abort("Error during initial size calculation (snprintf)\n");
    }
    // +1 for the null terminator
    r = snprintf((*com_prof)->name, requiredSize + 1, "%c%d.%d",(*com_prof)->id, id0, id1);
    if (r < 0) {
        mcpt_abort("Error during final formatting (snprintf)\n");
    } else if (r >= NAMELEN) {
            mcpt_abort("Name truncated during formatting\n");
    }
}

int
MPI_Pcontrol(const int level, ...)
{

    int mpi_errno = MPI_SUCCESS;
    va_list list;

    /* ... body of routine ...  */

    va_start(list, level);
    va_end(list);

    if ( level == 1)
        prof_enabled = 1;
    else if ( level == 0 )
        prof_enabled = 0;
    else
        printf("mpisee: MPI_Pcontrol called with invalid value: %d\nProfiling enabled = %d\n",level,prof_enabled);
    /* ... end of body of routine ... */
    return mpi_errno;
}

std::vector<std::string> convertToArrayOfStrings(char *proc_names, int size,
                                                 int name_length) {
    std::vector<std::string> machineNames;
    for (int i = 1; i < size; ++i) {
        std::string machineName(proc_names + i * name_length);
        machineNames.push_back(machineName);
    }
    return machineNames;
}

std::vector<std::string> convertToArrayOfPrims() {

    std::vector<std::string> machineNames;
    for (int i = 1; i < NUM_OF_PRIMS; ++i) {
        std::string PrimName(prim_names[i]);
        machineNames.push_back(PrimName);
    }
    return machineNames;
}


int
_MPI_Init(int *argc, char ***argv){
    int ret,rank,size;
    const auto start{std::chrono::steady_clock::now()};
    ret = PMPI_Init(argc, argv);
    const auto end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> duration{end - start};
    init_time = duration.count();
    PMPI_Comm_rank(MPI_COMM_WORLD, &rank);
    PMPI_Comm_size(MPI_COMM_WORLD, &size);


    MPI_Comm_create_keyval(MPI_COMM_DUP_FN,MPI_COMM_NULL_DELETE_FN,
                           &keys[0],NULL);

    MPI_Comm_create_keyval(MPI_COMM_DUP_FN,MPI_COMM_NULL_DELETE_FN,
                           &keys[1],NULL);

    if ( rank == 0 ){
        appname = (char*)malloc(sizeof(char)*1024);
        appname = get_appname();
        printf("MPI_Init: mpisee Profiling Tool version %d.%d\nProfiling application\
 %s\n",MPISEE_MAJOR_VERSION,MPISEE_MINOR_VERSION,appname);
 #ifdef MPICH_NAME
        printf("MPICH library used\n");
 #endif
 #ifdef OMPI_MAJOR_VERSION
        printf("OpenMPI library used\n");
 #endif
        fflush(stdout);
    }

    alloc_init_commprof(MPI_COMM_WORLD, 'W');
    comms_table.push_back(MPI_COMM_WORLD);
    profile_this(MPI_COMM_WORLD, 0, MPI_DATATYPE_NULL, Init, init_time, 0);

    if ( argc != NULL )
        ac = *argc;

    total_time = MPI_Wtime();

    return ret;
}


static int
_MPI_Init_thread(int *argc, char ***argv, int required, int *provided){
    int ret,rank,size;
    const auto start{std::chrono::steady_clock::now()};
    ret = PMPI_Init_thread(argc, argv, required, provided);
    const auto end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> duration{end - start};
    init_time = duration.count();
    PMPI_Comm_rank(MPI_COMM_WORLD, &rank);
    PMPI_Comm_size(MPI_COMM_WORLD, &size);


    MPI_Comm_create_keyval(MPI_COMM_DUP_FN,MPI_COMM_NULL_DELETE_FN,
                           &keys[0],NULL);

    MPI_Comm_create_keyval(MPI_COMM_DUP_FN,MPI_COMM_NULL_DELETE_FN,
                           &keys[1],NULL);

    if ( rank == 0 ){
        appname = (char*)malloc(sizeof(char)*1024);
        appname = get_appname();
        printf("MPI_Init_thread: mpisee Profiling Tool\nProfiling\
application %s\n",appname);
        fflush(stdout);
 #ifdef MPICH_NAME
        printf("MPICH library used\n");
 #endif
 #ifdef OMPI_MAJOR_VERSION
        printf("OpenMPI library used\n");
 #endif
        fflush(stdout);
    }

    alloc_init_commprof(MPI_COMM_WORLD, 'W');
    comms_table.push_back(MPI_COMM_WORLD);
    profile_this(MPI_COMM_WORLD, 0, MPI_DATATYPE_NULL, Init_thread, init_time, 0);


    if ( argc != NULL )
        ac = *argc;

    total_time = MPI_Wtime();

    return ret;
}


int
MPI_Init_thread(int *argc, char ***argv, int required, int *provided)
{
    if ( argc != NULL )
        getProcCmdLine (&ac, av);
    return _MPI_Init_thread(argc, argv, required, provided);
}


extern "C" {
void
F77_MPI_INIT_THREAD (int *required, int *provided, int *ierr)
{
    char **tmp;
    int ret;
    getProcCmdLine (&ac, av);
    tmp = av;
    ret = _MPI_Init_thread(&ac, (char***)&tmp , *required, provided);
    mpisee_fortran_in_place_init();
    *ierr = ret;
    return;
}
}

/* void
mpi_init_ (int *ierr){
  int ret = 0;
  char **tmp;
  getProcCmdLine (&ac, av);
  tmp = av;
  ret = _MPI_Init (&ac, (char ***) &tmp);
  *ierr = ret;
  return;
} */

extern "C" {
void
F77_MPI_INIT (int *ierr)
{
  int ret = 0;
  char **tmp;
  getProcCmdLine (&ac, av);
  tmp = av;
  ret = _MPI_Init (&ac, (char ***) &tmp);
  mpisee_fortran_in_place_init();
  *ierr = ret;
  return;
}
}

int
MPI_Init(int *argc, char ***argv)
{
    if ( argc != NULL  ){
        getProcCmdLine (&ac, av);
        //getRunCmd(ac, av);

    }
    return _MPI_Init(argc, argv);
}

/*
 * The naming scheme for the MPI_Comm_create is done as MPI_Comm_split
 * Uses the unique character 'c' as prefix
 */
int
MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm)
{
    int ret;
    ret = PMPI_Comm_create(comm, group, newcomm);
    if ( newcomm == NULL || *newcomm == MPI_COMM_NULL )
        return ret;

    alloc_init_commprof(*newcomm, 'c');
    comms_table.push_back(*newcomm);
    return ret;
}


extern "C" {
void
F77_MPI_COMM_CREATE(MPI_Fint  * comm, MPI_Fint  * group, MPI_Fint  *comm_out ,
                    MPI_Fint *ierr)
{
    int ret;
    MPI_Comm c_comm, c_comm_out;
    MPI_Group c_group;
    c_comm = MPI_Comm_f2c(*comm);
    c_group = MPI_Group_f2c(*group);

    ret= MPI_Comm_create(c_comm, c_group, &c_comm_out);
    *ierr = ret;
    if( ret == MPI_SUCCESS )
        *comm_out = MPI_Comm_c2f(c_comm_out);
    return;

}
}

/*
 * The naming scheme for the MPI_Comm_split requires two numbers
 * 1. current number of communicators
 * 2. a unique id for each new communicator via split
 * Uses the unique character 's' as prefix
 */
int
MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm)
{
    int ret;
    ret = PMPI_Comm_split(comm, color, key, newcomm);
    if ( newcomm == NULL || *newcomm == MPI_COMM_NULL )
        return ret;

    alloc_init_commprof(*newcomm, 's');
    comms_table.push_back(*newcomm);
    return ret;
}

extern "C" {
void
F77_MPI_COMM_SPLIT(MPI_Fint  * comm, int  * color, int  * key,
                   MPI_Fint  *comm_out , MPI_Fint *ierr)
{
    int ret;
    MPI_Comm c_comm,c_comm_out;
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Comm_split(c_comm, *color, *key, &c_comm_out);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *comm_out = MPI_Comm_c2f(c_comm_out);
    return;
}
}


int
MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm)
{
    int ret;
    ret = PMPI_Comm_dup(comm, newcomm);
    if ( newcomm == NULL || *newcomm == MPI_COMM_NULL )
        return ret;

    alloc_init_commprof(*newcomm, 'd');
    comms_table.push_back(*newcomm);
    return ret;
}


extern "C" {
void
F77_MPI_COMM_DUP(MPI_Fint  * comm, MPI_Fint  *comm_out , MPI_Fint *ierr)
{
    int ret;
    MPI_Comm newcomm;
    MPI_Comm c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Comm_dup(c_comm, &newcomm);
    *ierr = ret;
    if ( ret == MPI_SUCCESS  )
        *comm_out =  MPI_Comm_c2f(newcomm);
    return;
}
}

int
MPI_Comm_idup(MPI_Comm comm, MPI_Comm *newcomm, MPI_Request *request)
{

    int ret;
    ret = PMPI_Comm_idup(comm, newcomm, request);
    if ( newcomm == NULL || *newcomm == MPI_COMM_NULL )
        return ret;
    alloc_init_commprof(*newcomm, 'i');
    requests_map[*request] = comm;
    comms_table.push_back(*newcomm);
    return ret;
}

extern "C" {
void
mpi_comm_idup_(MPI_Fint  * comm, MPI_Fint  *comm_out, MPI_Fint  *request,
                 MPI_Fint *ierr)
{
    int ret;
    MPI_Comm c_comm, c_comm_out;
    MPI_Request c_request;
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Comm_idup(c_comm, &c_comm_out, &c_request);
    *ierr = ret;
    if ( ret == MPI_SUCCESS ){
        *comm_out = MPI_Comm_c2f(c_comm_out);
        *request = MPI_Request_c2f(c_request);
    }
    return;
}
}

int
MPI_Cart_create(MPI_Comm old_comm, int ndims, const int *dims,
                const int *periods, int reorder, MPI_Comm *newcomm)
{
    int ret;
    ret = PMPI_Cart_create(old_comm, ndims, dims, periods, reorder, newcomm);
    if ( newcomm == NULL || *newcomm == MPI_COMM_NULL ){
        return ret;
    }

    alloc_init_commprof(*newcomm, 'a');
    comms_table.push_back(*newcomm);
    return ret;
}


extern "C" {
void
F77_MPI_CART_CREATE(MPI_Fint  * comm_old, int  * ndims, const int  *dims,
                    const int  *periods, int  * reorder,
                    MPI_Fint  *comm_cart , MPI_Fint *ierr)
{
    int ret;
    MPI_Comm c_comm_old, c_comm_cart;
    c_comm_old = MPI_Comm_f2c(*comm_old);
    ret = MPI_Cart_create(c_comm_old, *ndims, dims, periods, *reorder, &c_comm_cart);
    *ierr = ret;
    if ( ret == MPI_SUCCESS  )
        *comm_cart = MPI_Comm_c2f(c_comm_cart);
    return;
}
}


int
MPI_Cart_sub(MPI_Comm comm, const int *remain_dims, MPI_Comm *newcomm)
{
    int ret;

    ret = PMPI_Cart_sub(comm, remain_dims, newcomm);
    if ( newcomm == NULL || *newcomm == MPI_COMM_NULL ){
        return ret;
    }
    alloc_init_commprof(*newcomm, 'b');
    comms_table.push_back(*newcomm);
    return ret;
}

extern "C" {
void
F77_MPI_CART_SUB(MPI_Fint  * comm, const int  *remain_dims,
                 MPI_Fint  *comm_new , MPI_Fint *ierr)
{
    int rc;
    MPI_Comm c_comm;
    MPI_Comm c_comm_new;
    c_comm = MPI_Comm_f2c(*comm);
    c_comm_new = MPI_Comm_f2c(*comm_new);
    rc = MPI_Cart_sub(c_comm, remain_dims, &c_comm_new);
    *comm_new = MPI_Comm_c2f(c_comm_new);
    *ierr = rc;
    return;
}
}

int
MPI_Graph_create(MPI_Comm comm_old, int nnodes, const int *index,
                 const int *edges, int reorder, MPI_Comm *newcomm)
{
    int ret;

    ret = PMPI_Graph_create(comm_old, nnodes, index, edges, reorder, newcomm);
    if ( newcomm == NULL || *newcomm == MPI_COMM_NULL ){
        return ret;
    }
    alloc_init_commprof(*newcomm, 'r');
    comms_table.push_back(*newcomm);
    return ret;
}


extern "C" {
void
F77_MPI_GRAPH_CREATE(MPI_Fint  * comm_old, int  * nnodes, const int  *index,
                     const int  *edges, int  * reorder, MPI_Fint  *newcomm,
                     MPI_Fint *ierr)
{
    int ret;
    MPI_Comm c_comm_old;
    MPI_Comm c_newcomm;

    c_comm_old = MPI_Comm_f2c(*comm_old);

    ret = MPI_Graph_create(c_comm_old, *nnodes, index, edges, *reorder, &c_newcomm);

    *ierr = (MPI_Fint)ret;
    if ( ret == MPI_SUCCESS ) {
        *newcomm = MPI_Comm_c2f(c_newcomm);
    }
    return;
}

}


int
MPI_Dist_graph_create(MPI_Comm comm_old, int n, const int *nodes,
                      const int *degrees, const int *targets,
                      const int *weights, MPI_Info info, int reorder,
                      MPI_Comm *newcomm)
{
    int ret;
    ret = PMPI_Dist_graph_create(comm_old, n, nodes, degrees, targets, weights, info, reorder, newcomm);
        if ( newcomm == NULL || *newcomm == MPI_COMM_NULL ){
        return ret;
    }
    alloc_init_commprof(*newcomm, 'g');
    comms_table.push_back(*newcomm);

    return ret;
}

extern "C" {
void
mpi_dist_graph_create_(MPI_Fint *comm_old, int *n, const int *nodes,
                               const int *degrees, const int *targets,
                               const int *weights, MPI_Fint *info, int *reorder,
                               MPI_Fint *newcomm, MPI_Fint *ierr)
{
    int ret;
    MPI_Comm c_comm_old;
    MPI_Comm c_newcomm;
    MPI_Info c_info;
    c_comm_old = MPI_Comm_f2c(*comm_old);
    c_info = MPI_Info_f2c(*info);
    ret = MPI_Dist_graph_create(c_comm_old, *n, nodes, degrees, targets, weights, c_info, *reorder, &c_newcomm);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *newcomm = MPI_Comm_c2f(c_newcomm);
    return;
}

}

int
MPI_Comm_split_type(MPI_Comm comm, int split_type, int key, MPI_Info info,
                    MPI_Comm *newcomm){
    int ret;
    ret = PMPI_Comm_split_type(comm, split_type, key, info, newcomm);
    if ( newcomm == NULL || *newcomm == MPI_COMM_NULL ){
        return ret;
    }
    alloc_init_commprof(*newcomm, 't');
    comms_table.push_back(*newcomm);
    return ret;
}

extern "C" {
void
F77_MPI_COMM_SPLIT_TYPE(MPI_Fint  * comm, int  * split_type, int  * key,
                        MPI_Fint *info, MPI_Fint  *newcomm , MPI_Fint *ierr)
{
    int ret;
    MPI_Comm c_comm,c_comm_out;
    MPI_Info c_info;
    c_info = PMPI_Info_f2c(*info);
    c_comm = PMPI_Comm_f2c(*comm);
    ret = MPI_Comm_split_type(c_comm, *split_type, *key, c_info,&c_comm_out);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *newcomm = MPI_Comm_c2f(c_comm_out);
    return;
}
}

int
MPI_Comm_create_group(MPI_Comm comm, MPI_Group group, int tag, MPI_Comm *newcomm)
{
    int ret;
    ret = PMPI_Comm_create_group(comm, group, tag, newcomm);
    if ( newcomm == NULL || *newcomm == MPI_COMM_NULL ){
        return ret;
    }
    alloc_init_commprof(*newcomm, 'u');
    comms_table.push_back(*newcomm);
    return ret;
}

extern "C" {
void
mpi_comm_create_group_(MPI_Fint *comm, MPI_Fint *group, int *tag,
                          MPI_Fint *comm_out , MPI_Fint *ierr)
{
    int ret;
    MPI_Comm c_comm, c_comm_out;
    MPI_Group c_group;
    c_comm = MPI_Comm_f2c(*comm);
    c_group = MPI_Group_f2c(*group);
    ret = MPI_Comm_create_group(c_comm, c_group, *tag, &c_comm_out);
    *ierr = ret;
    if( ret == MPI_SUCCESS )
        *comm_out = MPI_Comm_c2f(c_comm_out);
    return;
}
}

int
MPI_Dist_graph_create_adjacent(MPI_Comm comm_old, int indegree,
                               const int *sources, const int *sourceweights,
                               int outdegree, const int *destinations,
                               const int *destweights, MPI_Info info,
                               int reorder, MPI_Comm *comm_dist_graph)
{
    int ret;
    ret = PMPI_Dist_graph_create_adjacent(comm_old, indegree, sources,
                                          sourceweights, outdegree,
                                          destinations, destweights, info,
                                          reorder, comm_dist_graph);
    if ( comm_dist_graph == NULL || *comm_dist_graph == MPI_COMM_NULL ){
        return ret;
    }
    alloc_init_commprof(*comm_dist_graph, 'j');
    comms_table.push_back(*comm_dist_graph);
    return ret;
}

int
MPI_Wait(MPI_Request *request, MPI_Status *status)
{
    int ret;
    double t_elapsed;
    MPI_Comm comm;
    if ( prof_enabled == 1 ){
        auto it = requests_map.find(*request);
        if (it != requests_map.end()) {
            comm = requests_map[*request];
        }
        else {
            comm = MPI_COMM_WORLD;
        }

        t_elapsed = MPI_Wtime();
        ret = PMPI_Wait(request, status);
        t_elapsed = MPI_Wtime() - t_elapsed;
        if ( comm != MPI_COMM_NULL  ){
            profile_this(comm, 0, MPI_DATATYPE_NULL, Wait, t_elapsed, 1);
            requests_map.erase(*request);
        }
        else {
            mcpt_abort("mpisee: NULL COMMUNICATOR in MPI_Wait\n");
            return ret;
        }
    }
    else{
        ret = PMPI_Wait(request, status);
    }
    return ret;
}

extern "C" {
void
F77_MPI_WAIT(MPI_Fint  *request, MPI_Status  *status , MPI_Fint *ierr)
{
   int ret;
   MPI_Request c_request;
   c_request = MPI_Request_f2c(*request);
   ret = MPI_Wait(&c_request, status);
   *ierr = ret;
}
}


int
MPI_Waitall(int count, MPI_Request array_of_requests[],
            MPI_Status array_of_statuses[])
{
    int ret,i;
    double t_elapsed;
    MPI_Comm comm ;
    if ( prof_enabled == 1 ){
        for (i = 0; i < count; i++) {
            auto it = requests_map.find(array_of_requests[i]);
            if (it != requests_map.end()) {
                comm = requests_map[array_of_requests[i]];
                break;
            }
            else {
                comm = MPI_COMM_WORLD;
            }
        }
        t_elapsed = MPI_Wtime();
        ret = PMPI_Waitall(count, array_of_requests, array_of_statuses);
        t_elapsed = MPI_Wtime() - t_elapsed;
        if ( comm != MPI_COMM_NULL){
            profile_this(comm, 0, MPI_DATATYPE_NULL, Waitall, t_elapsed, 1);
            for (i = 0; i < count; i++) {
                requests_map.erase(array_of_requests[i]);
            }
        }
        else{
             mcpt_abort("NULL COMMUNICATOR in MPI_Waitall\n");
             return ret;
        }
    }
    else{
        ret = PMPI_Waitall(count, array_of_requests, array_of_statuses);
    }
    return ret;
}


extern "C" {
void
F77_MPI_WAITALL(int  * count, MPI_Fint  *array_of_requests,
                MPI_Status  *array_of_statuses , MPI_Fint *ierr)
{

   int ret,i;
   MPI_Request *c_requests;
   c_requests = (MPI_Request*) malloc (sizeof(MPI_Request)*(*count));
   for ( i =0; i<*count; i++ ){
       c_requests[i] = MPI_Request_f2c(array_of_requests[i]);
   }
   ret = MPI_Waitall(*count, c_requests, array_of_statuses);
   *ierr = ret;
   if ( ret == MPI_SUCCESS ){
       for ( i =0; i<*count; i++ ){
           array_of_requests[i] = MPI_Request_c2f(c_requests[i]);
       }
   }
   free( c_requests );
}
}

int
MPI_Waitany(int count, MPI_Request *array_of_requests, int *index, MPI_Status *status)
{

    int ret,i;
    double t_elapsed;
    MPI_Comm *comm_array;

    if ( prof_enabled == 1 ){
        comm_array = (MPI_Comm*) malloc (sizeof(MPI_Comm)*count);
        for ( i =0; i<count; i++ ){
            auto it = requests_map.find(array_of_requests[i]);
            if (it != requests_map.end()) {
                comm_array[i] = requests_map[array_of_requests[i]];
            }
            else {
                comm_array[i] = MPI_COMM_NULL;
            }
        }
        t_elapsed = MPI_Wtime();
        ret = PMPI_Waitany(count, array_of_requests,index,status);
        t_elapsed = MPI_Wtime() - t_elapsed;
        if ( comm_array[*index] != MPI_COMM_NULL){
            profile_this(comm_array[*index], 0, MPI_DATATYPE_NULL, Waitany, t_elapsed, 1);
            requests_map.erase(array_of_requests[*index]);

        }
        // else{
        //     mcpt_abort("mpisee: NULL COMMUNICATOR in MPI_Waitany\n");
        //     return ret;
        // }

        free(comm_array);
    }
    else{
        ret = PMPI_Waitany(count, array_of_requests,index,status);
    }
    return ret;
}

extern "C" {
void
F77_MPI_WAITANY(int  * count, MPI_Fint  *array_of_requests, int  *index,
                MPI_Status  *status , MPI_Fint *ierr)
{
    int ret,i;
    MPI_Request *c_array_of_requests = NULL;
    //c_array_of_requests = (MPI_Request*)malloc(sizeof(MPI_Request)*(*count));
    //assert(c_array_of_requests);
    for (i = 0; i < *count; i++) {
        c_array_of_requests[i] = MPI_Request_f2c(array_of_requests[i]);
    }

    ret = MPI_Waitany(*count, c_array_of_requests, index, status);

    *ierr = (MPI_Fint)ret;
    free(c_array_of_requests);
    return;

}
}

int
MPI_Test(MPI_Request *request, int *flag, MPI_Status *status)
{

    int ret;
    double t_elapsed;
    MPI_Comm comm;

    if ( prof_enabled == 1 ){
        auto it = requests_map.find(*request);
        if (it != requests_map.end()) {
            comm = requests_map[*request];
        }
        else {
            comm = MPI_COMM_NULL;
        }
        t_elapsed = MPI_Wtime();
        ret = PMPI_Test(request,flag,status);
        t_elapsed = MPI_Wtime() - t_elapsed;
        if ( comm != MPI_COMM_NULL ){
            profile_this(comm, 0, MPI_DATATYPE_NULL, Test, t_elapsed, 1);
            if ( *flag == 1 ){
                requests_map.erase(*request);
            }
        }
        else{
            return ret;
        }
    }
    else{
        ret = PMPI_Test(request,flag,status);
    }
    return ret;
}


extern "C" {
void
F77_MPI_TEST(MPI_Fint  *request, int  *flag, MPI_Status  *status , MPI_Fint *ierr)
{
    int ret;
    MPI_Request c_request;

    c_request = MPI_Request_f2c(*request);
    ret = MPI_Test(&c_request,flag,status);
    *ierr = (MPI_Fint)ret;
    if ( ret == MPI_SUCCESS ) {
        *request = MPI_Request_c2f(c_request);
    }
    return;
}
}

int
MPI_Testany(int count, MPI_Request *array_of_requests, int *index, int *flag, MPI_Status *status)
{


    int ret,i;
    double t_elapsed;
    MPI_Comm *comm_array;

    if ( prof_enabled == 1 ){
        comm_array = (MPI_Comm*) malloc (sizeof(MPI_Comm)*count);
        for ( i =0; i<count; i++ ){
            auto it = requests_map.find(array_of_requests[i]);
            if (it != requests_map.end()) {
                comm_array[i] = requests_map[array_of_requests[i]];
            }
            else {
                comm_array[i] = MPI_COMM_NULL;
            }
        }
        t_elapsed = MPI_Wtime();
        ret = PMPI_Waitany(count, array_of_requests,index,status);
        t_elapsed = MPI_Wtime() - t_elapsed;
        if ( comm_array[*index] != MPI_COMM_NULL){
            profile_this(comm_array[*index], 0, MPI_DATATYPE_NULL, Testany, t_elapsed, 1);
            if ( *flag == 1 ){
                requests_map.erase(array_of_requests[*index]);
            }
        }

        free(comm_array);
    }
    else{
        ret = PMPI_Testany(count, array_of_requests, index, flag, status);
    }
    return ret;

}


extern "C" {
void
F77_MPI_TESTANY(int  * count, MPI_Fint  *array_of_requests, int  *index,
                int  *flag, MPI_Status  *status , MPI_Fint *ierr)
{

    int ret,i;
    MPI_Request *c_array_of_requests = NULL;
    //c_array_of_requests = (MPI_Request*)malloc(sizeof(MPI_Request)*(*count));
    //assert(c_array_of_requests != NULL);
    for (i = 0; i < *count; i++) {
        c_array_of_requests[i] = MPI_Request_f2c(array_of_requests[i]);
    }

    ret = MPI_Testany(*count, c_array_of_requests, index, flag, status);

    *ierr = (MPI_Fint)ret;
    free(c_array_of_requests);
    return;

}
}



int
MPI_Comm_free(MPI_Comm *comm)
{
    int ret,flag;
    prof_metadata *metadata;
    comm_profiler *data;
    MPI_Group group;

    PMPI_Comm_get_attr(*comm, keys[0], &metadata, &flag);
    // Debug prints with flag check
    /*
    if (flag) {
        std::cout << "mpisee: Comm_free: Comm_get_attr metadata found in communicator\n";
    } else {
        mcpt_abort("Comm_free: Comm_get_attr did not find metadata in communicator\n");
    }
     */

    PMPI_Comm_get_attr(*comm, keys[1], &data, &flag);
    // Debug prints with flag check
    /*
    if (flag) {
        std::cout << "mpisee: Comm_free: Comm_get_attr data found in data\n";
    } else {
        mcpt_abort("Comm_free: Comm_get_attr d data in communicator\n");
    }
     */

    prof_meta_pair *free_pair = new prof_meta_pair();
    free_pair->meta = *metadata;
    free_pair->prof = *data;
    MPI_Comm_group(*comm, &group);

    free_array.push_back(std::make_pair(free_pair, group));

    // Find newcomm in comms and remove it
    auto it = std::find(comms_table.begin(), comms_table.end(), *comm);
    if (it != comms_table.end()) {
        comms_table.erase(it);
    } else {
        mcpt_abort("Comm_free: Comm not found in comms_table\n");
    }
    comms_table.shrink_to_fit();

    ret = PMPI_Comm_free(comm);
    return ret;
}


extern "C" {
void
F77_MPI_COMM_FREE(MPI_Fint *comm, MPI_Fint *ierr)
{
    int ret;
    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Comm_free(&c_comm);
    *ierr = ret;
    return;
}
}

static int
_Finalize(void) {
    int rank, size, buf[2] = {-1, -1};
    int flag, len, resultlen, datasize;
    std::vector<comm_all> recv_comm_buffer;
    long unsigned num_of_comms;
    char version[MPI_MAX_LIBRARY_VERSION_STRING];
    char proc_name[MPI_MAX_PROCESSOR_NAME];
    char *proc_names = NULL;
    double *alltimes = NULL;
    std::vector<double> mpi_times;
    long unsigned total_num_of_comms = 0;
    std::vector<prof_attrs*> local_communicators;
    // Place the profiling data into a vector and gather it to rank 0
    std::vector<comm_all> metadata_array;
    std::vector<comm_data> data_array;
    comm_all comm_meta;
    prof_metadata *metadata;
    comm_profiler *comm_prof_data;
    MPI_Comm newcomm;

    total_time = MPI_Wtime() - total_time;


    if (PMPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
        mcpt_abort("MPI_Barrier failed\n");
    }


    PMPI_Comm_rank(MPI_COMM_WORLD, &rank);
    PMPI_Comm_size(MPI_COMM_WORLD, &size);


    // Re-create the communicators from the group_table
    if ( rank == 0)
        std::cout << "mpisee: group_table size = " << free_array.size() << std::endl;
    for (long unsigned i = 0; i < free_array.size(); ++i) {
        PMPI_Comm_create_group(MPI_COMM_WORLD, free_array[i].second, 0, &newcomm);
        if (newcomm == MPI_COMM_NULL) {
            mcpt_abort("Comm_create_group failed\n");
        }
        PMPI_Comm_set_attr(newcomm, keys[0], &free_array[i].first->meta);
        PMPI_Comm_set_attr(newcomm, keys[1], &free_array[i].first->prof);
        comms_table.push_back(newcomm);
    }
    num_of_comms = comms_table.size();

    for(long unsigned i = 0; i < num_of_comms; ++i) {
        PMPI_Comm_get_attr(comms_table[i], keys[0], &metadata, &flag);
        buf[0] = rank;
        buf[1] = metadata->comms;
        PMPI_Bcast(buf, 2, MPI_INT, 0, comms_table[i]);
        overwrite_name(&metadata, buf[0], buf[1]);
        strcpy(comm_meta.name, metadata->name);
        comm_meta.size = metadata->size;
        datasize = 0;
        PMPI_Comm_get_attr(comms_table[i], keys[1], &comm_prof_data, &flag);
        if (!flag) {
            mcpt_abort("Finalize: Comm_get_attr failed\n");
        }
        for (int k=0; k<NUM_OF_PRIMS; ++k) {
            for (int j = 0; j < NUM_BUCKETS; ++j) {
                // Check if the key exists in the map

                //if (rank == 0)
                //    std::cout << "mpisee: key = " << key << std::endl;
                auto it = comm_prof_data->map.find(getPrimBucketKey(k, j));

                if (it != comm_prof_data->map.end()) {
                    // Allocate comm_data struct and copy the data from map
                    comm_data data;
                    data.comm_id = i;
                    data.prim = k;
                    data.bucketIndex = j;
                    data.num_messages = it->second.num_messages;
                    data.time = it->second.time;
                    data.volume = it->second.volume;
                    data_array.push_back(data);
                    datasize++;
                }
            }
        }
        comm_meta.datasize = datasize;
        //std::cout << "mpisee: comm_meta.datasize = " << comm_meta.datasize << std::endl;
        metadata_array.push_back(comm_meta); // metadata array send seperately
    }


    int *c_recvcounts = NULL;
    int *c_displs = NULL;

    if ( rank == 0 ){
        c_recvcounts = (int *)malloc(sizeof(int) * size);
        if (c_recvcounts == NULL) {
            mcpt_abort("malloc error for recvcounts");
        }
        c_displs = (int *)malloc(sizeof(int) * size);
        if (c_displs == NULL) {
            mcpt_abort("malloc error for displs");
        }

    }
    PMPI_Gather(&num_of_comms, 1, MPI_INT, c_recvcounts,
                1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // Compute displacements
        c_displs[0] = 0;
        total_num_of_comms = c_recvcounts[0];
        for (int i = 1; i < size; ++i) {
          c_displs[i] = c_displs[i - 1] + c_recvcounts[i - 1];
          total_num_of_comms += c_recvcounts[i];
        }
        std::cout << "mpisee: total number of communicators = " << total_num_of_comms << std::endl;
        recv_comm_buffer.resize(total_num_of_comms);

    }

    // Gather the metadata
    comm_all dummy;
    MPI_Datatype MPI_COMM_ALL;
    MPI_Datatype types[3] = {MPI_CHAR, MPI_INT, MPI_INT};
    int blocklengths[3] = {NAMELEN, 1, 1};
    MPI_Aint displacements[3];
    MPI_Aint base;
    PMPI_Get_address(&dummy, &base);
    PMPI_Get_address(&dummy.name, &displacements[0]);
    PMPI_Get_address(&dummy.size, &displacements[1]);
    PMPI_Get_address(&dummy.datasize, &displacements[2]);
    // Convert addresses to displacements
    for (int i = 0; i < 3; i++) {
        displacements[i] = MPI_Aint_diff(displacements[i], base);
    }

    PMPI_Type_create_struct(3, blocklengths, displacements, types, &MPI_COMM_ALL);
    PMPI_Type_commit(&MPI_COMM_ALL);

    PMPI_Gatherv(metadata_array.data(), num_of_comms, MPI_COMM_ALL,
                recv_comm_buffer.data(), c_recvcounts, c_displs, MPI_COMM_ALL, 0, MPI_COMM_WORLD);


    PMPI_Type_free(&MPI_COMM_ALL);
    metadata_array.clear();
    metadata_array.shrink_to_fit(); //recv_comm_buffer has all communicator metadata now

    // Gather the actual data now
    // 1. Create MPI datatype for comm_data
    comm_data dummy_data;
    MPI_Datatype MPI_COMM_DATA;
    MPI_Datatype datatypes[6] = {MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_UINT64_T};
    int blocklengths2[6] = {1, 1, 1, 1, 1, 1};
    MPI_Aint displacements2[6]; // Use more descriptive names for clarity
    MPI_Aint comm_data_base;    // More descriptive name

    PMPI_Get_address(&dummy_data, &comm_data_base);

    PMPI_Get_address(&dummy_data.comm_id, &displacements2[0]);
    PMPI_Get_address(&dummy_data.prim, &displacements2[1]);
    PMPI_Get_address(&dummy_data.bucketIndex, &displacements2[2]);
    PMPI_Get_address(&dummy_data.time, &displacements2[3]);
    PMPI_Get_address(&dummy_data.num_messages, &displacements2[4]);
    PMPI_Get_address(&dummy_data.volume, &displacements2[5]);

    for (int i = 0; i < 6; ++i) {
        displacements2[i] = MPI_Aint_diff(displacements2[i], comm_data_base);
    }

    PMPI_Type_create_struct(6, blocklengths2, displacements2, datatypes, &MPI_COMM_DATA);

    PMPI_Type_commit(&MPI_COMM_DATA);

    // 2. Gather the profiling data from all ranks to rank 0
    int local_data_size = data_array.size();
    //std::cout << "mpisee: local_data_size = " << local_data_size << std::endl;
    int total_num_of_data = 0;
    int *recvcounts = NULL;
    int *displs = NULL;
    // allocate recvcounts and displs only for rank 0
    if (rank == 0) {
        recvcounts = (int *)malloc(sizeof(int) * size);
        if (recvcounts == NULL) {
            mcpt_abort("malloc error for recvcounts buffer Rank: %d\n", rank);
        }
        displs = (int *)malloc(sizeof(int) * size);
        if (displs == NULL) {
            mcpt_abort("malloc error for displs buffer Rank: %d\n", rank);
        }
    }
    PMPI_Gather(&local_data_size, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Root process: Allocate receiving buffer and calculate displacements
    std::vector<comm_data> recv_data_buffer;
    if (rank == 0) {
        displs[0] = 0;
        total_num_of_data = recvcounts[0];
        for (int i = 1; i < size; ++i) {
            displs[i] = displs[i - 1] + recvcounts[i - 1];
            total_num_of_data += recvcounts[i];
        }
        //std::cout << "mpisee: total number of data = " << total_num_of_data << std::endl;
        recv_data_buffer.resize(total_num_of_data);
    }
    // Gather the data
    PMPI_Gatherv(data_array.data(), local_data_size, MPI_COMM_DATA,
                recv_data_buffer.data(), recvcounts, displs, MPI_COMM_DATA, 0, MPI_COMM_WORLD);

    // Clear the data_array vector
    data_array.clear();
    data_array.shrink_to_fit(); // recv_data_buffer has all data now


    PMPI_Get_processor_name(proc_name, &len);

    if (rank == 0) {
        proc_names = (char *)malloc(sizeof(char) * MPI_MAX_PROCESSOR_NAME * size);
        if (proc_names == NULL) {
          mcpt_abort("malloc error for proc_names buffer Rank: %d\n", rank);
        }
        alltimes = (double *)malloc(sizeof(double) * size);
        if (alltimes == NULL) {
          mcpt_abort("malloc error for alltimes buffer Rank: %d\n", rank);
        }
     }

    PMPI_Gather(proc_name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, proc_names,
                 MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 0, MPI_COMM_WORLD);


    PMPI_Gather(&total_time, 1, MPI_DOUBLE, alltimes, 1, MPI_DOUBLE, 0,
                MPI_COMM_WORLD);

    if ( rank == 0 ){
        int rc,commId,maxsize,minsize;
        sqlite3 *db = NULL;
        std::string outfile;
        int proc;
        double t;
        const char *env_var = getenv("MPISEE_OUTFILE");
        if (env_var != NULL) {
          rc = sqlite3_open_v2(env_var, &db, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_FULLMUTEX, NULL);
          if (rc) {
              mcpt_abort("Can't open database: error: %s\n",env_var, sqlite3_errmsg(db));
          } else {
              std::cout << "mpisee: Opened database: " << env_var << " successfully " << std::endl;
          }
          outfile = env_var;
        }
        else{
            int maxRetries = 3600;
            db = openSQLiteDBExclusively("mpisee_", ".db", maxRetries, outfile);
            if (db == NULL) {
                mcpt_abort("Error: Failed to open SQLite database exclusively after %d retries", maxRetries);
                return 1;
            }
        }

        PMPI_Get_library_version(version, &resultlen);
        for (size_t i = 0; i < strlen(version); i++)
        {
            if(version[i] == '\n')
            {
                version[i] = ' ';
            }
        }

        t = MPI_Wtime();
        createTables(db);
        std::cout << "mpisee: Writing the metadata table" << std::endl;

        insertMetadata(db, version, size, av, ac, MPISEE_MAJOR_VERSION,
                       MPISEE_MINOR_VERSION, mpisee_build_date,
                       mpisee_build_time, env_var);

        std::cout << "mpisee: Writing the MPI operations table" << std::endl;


        insertIntoOperationsEmpty(db, prim_names[0]);
        std::vector<std::string> operations = convertToArrayOfPrims();
        BatchInsertIntoOperations(db, operations);
        operations.clear();
        operations.shrink_to_fit();

        std::vector<double> times;
        if (alltimes != NULL){
          std::cout << "mpisee: Writing the exectimes table" << std::endl;
          insertIntoTimes(db, alltimes[0]);
          for (int i = 1; i < size; i++) {
            times.push_back(alltimes[i]);
          }
          BatchInsertIntoTimes(db, times);
          times.clear();
          times.shrink_to_fit();
          free(alltimes);
        } else {
          std::cout << "mpisee: Execution times NULL" << std::endl;
        }

        std::string machineName(proc_name);
        insertIntoMappings(db, machineName);
        std::vector<std::string> machines =
            convertToArrayOfStrings(proc_names, size, MPI_MAX_PROCESSOR_NAME);
        BatchInsertIntoMappings(db, machines);
        machines.clear();
        machines.shrink_to_fit();
        free(proc_names);


        std::vector<CommData> comms;
        std::vector<int> commIds;

        for (long unsigned i = 0; i < total_num_of_comms; i++) {
            // debug print
            //std::cout << "mpisee: Writing metadata for communicator: "
            // << recv_comm_buffer[i].name << ", size: " << recv_comm_buffer[i].size
            // << ", datasize: " << recv_comm_buffer[i].datasize << std::endl;
            comms.push_back({recv_comm_buffer[i].name, recv_comm_buffer[i].size});
        }

        commIds=CommsInsert(db, comms);
        if (commIds.size() != total_num_of_comms) {
            mcpt_abort("mpisee: Error: CommIds size does not match total_num_of_comms\n");
        }
        comms.clear();
        comms.shrink_to_fit();

        std::vector<DataEntry> entries;
        std::cout << "mpisee: Writing the main data table"
                  << std::endl;
        commId = 0;
        int comms_per_proc;
        int datalen,index = 0;
        for (proc = 0; proc < size; ++proc) {
            comms_per_proc = c_recvcounts[proc];
            for (int i = 0; i < comms_per_proc; ++i) {
                commId = commIds[c_displs[proc]+i];
                datalen = recv_comm_buffer[c_displs[proc]+i].datasize;
                for (int j = 0; j < datalen; ++j) {
                    //std::cout << "mpisee: Writing data for communicator: "
                    //          << commId << ", datasize: " << datalen
                    //          << ", index: " << index << std::endl;

                    if (recv_data_buffer[index].bucketIndex == 0) {
                        minsize = 0;
                        maxsize = buckets[0];
                    } else if (recv_data_buffer[index].bucketIndex == NUM_BUCKETS - 1) {
                        minsize = buckets[NUM_BUCKETS - 2];
                        maxsize = INT_MAX;
                    } else {
                        minsize = buckets[recv_data_buffer[index].bucketIndex - 1];
                        maxsize = buckets[recv_data_buffer[index].bucketIndex];
                    }
                    // Debug print
                    //std::cout << "mpisee: Writing data for communicator: " << commId << ", prim: "
                    // << recv_data_buffer[i].prim << ", minsize: " << minsize
                    // << ", maxsize: " << maxsize << ", num_messages: "
                    // << recv_data_buffer[i].num_messages << ", time: "
                    // << recv_data_buffer[i].time << ", volume: "
                    // << recv_data_buffer[i].volume << std::endl;
                    insertIntoDataEntry(entries, proc, commId, recv_data_buffer[index].prim,
                                        minsize, maxsize, recv_data_buffer[index].num_messages,
                                        recv_data_buffer[index].time, recv_data_buffer[index].volume);
                    index++;
                }
            }
        }
        executeBatchInsert(db, entries);
        t = MPI_Wtime() - t;
        std::cout << "mpisee: Output database file: " << outfile << ", time to write: " << t << " seconds" << std::endl;
        sqlite3_close(db);
        // Free buffers allocated by rank 0 only
//        free(c_displs);
//        free(displs);
//        free(c_recvcounts);
//        free(recvcounts);
//        free(proc_names);
//        free(alltimes);
    }

    PMPI_Barrier(MPI_COMM_WORLD);
    PMPI_Type_free(&MPI_COMM_DATA);
    return PMPI_Finalize();
}


int
MPI_Finalize (void)
{
  int rc = 0;

  rc = _Finalize ();

  return rc;
}

extern "C" {
void
F77_MPI_FINALIZE (int *ierr)
{
  int rc = 0;

  rc = _Finalize ();
  *ierr = rc;

  return;
}
}

