#ifndef COMMPROF_H_
#define COMMPROF_H_

#define MAX_ARGS 1024
#define MAX_DIMS 8
#include "utils.h"
#include "datastructlib/table.h"

extern int prof_enabled;
extern prof_attrs **local_data;
extern prof_attrs **local_comms;
extern int local_cid;
extern int my_coms;
extern Table_T request_tab;

int namekey();
int namedel(MPI_Comm comm, int keyval, void *attr, void *s);

#endif
