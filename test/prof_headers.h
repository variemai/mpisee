#ifndef PROF_HEADERS_H_
#define PROF_HEADERS_H_

#include <string.h>
extern int line_called;
extern char file_called[32];

#define MPI_Comm_split(comm,color,key,newcomm)                          \
    do { line_called = __LINE__;                                        \
        strcpy(file_called, __FILE__);                                  \
        MPI_Comm_split(comm,color,key,newcomm);  }                       \
    while (0)

#define MPI_Comm_create(comm,group,newcomm)                             \
    do { line_called = __LINE__;                                        \
        strcpy(file_called, __FILE__);                                  \
        MPI_Comm_create(comm,group,newcomm);  }                          \
    while (0)
#endif // PROF_HEADERS_H_
