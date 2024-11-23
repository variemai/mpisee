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
#include "include/commprof.h"
#include <mpi.h>
#include "include/mpisee_fortran.h"

int
MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root,
          MPI_Comm comm)
{
    int ret,rank;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Bcast(buffer, count, datatype, root, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;
        PMPI_Comm_rank(comm, &rank);
        profile_this(comm,count,datatype,Bcast,t_elapsed,0);
    }
    else{
        ret = PMPI_Bcast(buffer, count, datatype, root, comm);
    }
    return ret;

}

extern "C" {
void
mpi_bcast_(void  *buffer, int  * count, MPI_Fint  * datatype, int  * root,
              MPI_Fint  * comm , MPI_Fint *ierr)
{

    int ret;
    MPI_Datatype c_datatype;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);

    if ( buffer == mpisee_fortran_mpi_bottom )
        buffer = MPI_BOTTOM;

    ret = MPI_Bcast(buffer, *count, c_datatype, *root, c_comm);

    *ierr = ret;
}
}

int
MPI_Ibcast(void *buffer, int count, MPI_Datatype datatype, int root,
           MPI_Comm comm, MPI_Request *request)
{
    int ret,rank;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Ibcast(buffer, count, datatype, root, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;

        PMPI_Comm_rank(comm, &rank);
        profile_this(comm,count,datatype,Ibcast,t_elapsed,0);
        requests_map[*request] = comm;
    }
    else{

        ret = PMPI_Ibcast(buffer, count, datatype, root, comm, request);
    }
    return ret;
}

extern "C" {
void
mpi_ibcast_(void  *buffer, int  * count, MPI_Fint  * datatype, int  * root,
               MPI_Fint  * comm, MPI_Fint  *request , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    MPI_Request c_request;
    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);

    if ( buffer == mpisee_fortran_mpi_bottom )
        buffer = MPI_BOTTOM;

    ret = MPI_Ibcast(buffer, *count, c_datatype, *root, c_comm,&c_request);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
}
}

int
MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
    int ret;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;

        profile_this(comm,count,datatype,Allreduce,t_elapsed,0);
    }
    else{
        ret = PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
    }
    return ret;
}

extern "C" {
void
mpi_allreduce_(void  *sendbuf, void  *recvbuf, int  * count,
                  MPI_Fint  * datatype, MPI_Fint  * op, MPI_Fint  * comm ,
                  MPI_Fint *ierr)
{
    int ret;
    MPI_Op c_op;
    MPI_Comm c_comm;
    MPI_Datatype c_datatype;

    c_datatype = MPI_Type_f2c(*datatype);
    c_op = MPI_Op_f2c(*op);
    c_comm = MPI_Comm_f2c(*comm);

    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Allreduce(sendbuf, recvbuf, *count, c_datatype, c_op, c_comm);

    *ierr =ret;
    return;
}
}

int
MPI_Iallreduce(const void *sendbuf, void *recvbuf, int count,
               MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
               MPI_Request *request)
{
    int ret;
    double t_elapsed;
    if (prof_enabled == 1){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Iallreduce(sendbuf, recvbuf, count, datatype, op, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;

        profile_this(comm, count, datatype, Iallreduce, t_elapsed, 0);
        requests_map[*request] = comm;
    }
    else{
        ret = PMPI_Iallreduce(sendbuf, recvbuf, count, datatype, op, comm, request);
    }
    return ret;
}


extern "C" {
void
mpi_iallreduce_(void  *sendbuf, void  *recvbuf, int  * count,
                   MPI_Fint  * datatype, MPI_Fint  * op, MPI_Fint  * comm,
                   MPI_Fint  *request , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Op c_op;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_datatype = MPI_Type_f2c(*datatype);
    c_op = MPI_Op_f2c(*op);
    c_comm = MPI_Comm_f2c(*comm);


    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Iallreduce(sendbuf, recvbuf, *count, c_datatype, c_op, c_comm, &c_request);


    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
    return;
}
}

int
MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
              void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
{
    int ret;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                             recvtype, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;

        profile_this(comm,recvcount,recvtype,Allgather,t_elapsed,0);
    }
    else{
        ret = PMPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                             recvtype, comm);
    }
    return ret;
}


extern "C" {
void
mpi_allgather_(void *sendbuf, int  * sendcount, MPI_Fint  * sendtype,
                  void  *recvbuf, int  * recvcount, MPI_Fint  * recvtype,
                  MPI_Fint  * comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Comm c_comm;
    MPI_Datatype c_sendtype, c_recvtype;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;


    ret = MPI_Allgather(sendbuf, *sendcount, c_sendtype, recvbuf, *recvcount, c_recvtype, c_comm);

    *ierr =ret;
    return;
}
}

int
MPI_Iallgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
              void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
               MPI_Request *request)
{
    int ret;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Iallgather(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                             recvtype, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;

        profile_this(comm,recvcount,recvtype,Iallgather,t_elapsed,0);
        requests_map[*request] = comm;
    }
    else{
        ret = PMPI_Iallgather(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                              recvtype, comm, request);
    }
    return ret;
}

extern "C" {
void
mpi_iallgather_(void *sendbuf, int  * sendcount, MPI_Fint  * sendtype,
                           void  *recvbuf, int  * recvcount, MPI_Fint  * recvtype,
                           MPI_Fint  * comm , MPI_Fint  *request , MPI_Fint *ierr)
{
    int ret;
    MPI_Comm c_comm;
    MPI_Datatype c_sendtype, c_recvtype;
    MPI_Request c_request;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);


    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Iallgather(sendbuf, *sendcount, c_sendtype, recvbuf, *recvcount, c_recvtype, c_comm, &c_request);

    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
    return;
}

}


int
MPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
             void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
{
    int ret;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                            recvtype, comm);

        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm,sendcount,sendtype,Alltoall,t_elapsed,0);
    }
    else{
        ret = PMPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                            recvtype, comm);
    }
    return ret;
}

extern "C" {
void
mpi_alltoall_(void  *sendbuf, int  * sendcount, MPI_Fint  * sendtype,
                 void  *recvbuf, int  * recvcnt, MPI_Fint  * recvtype,
                 MPI_Fint  * comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    // All processes must use MPI_IN_PLACE
    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Alltoall(sendbuf, *sendcount, c_sendtype, recvbuf, *recvcnt, c_recvtype, c_comm);

    *ierr = ret;
    return;
}
}

int
MPI_Ialltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                  void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                  MPI_Request *request)
{
    int ret;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Ialltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                             recvtype, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm,sendcount,sendtype,Ialltoall,t_elapsed,0);
        requests_map[*request] = comm;
    }
    else{
        ret = PMPI_Ialltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                             recvtype, comm, request);
    }
    return ret;
}

extern "C" {
void
mpi_ialltoall_(void *sendbuf, int *sendcount, MPI_Fint *sendtype,
                       void *recvbuf, int *recvcount, MPI_Fint  *recvtype,
                       MPI_Fint  *comm , MPI_Fint *request , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

   // All processes must use MPI_IN_PLACE
    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Ialltoall(sendbuf, *sendcount, c_sendtype, recvbuf,
                        *recvcount, c_recvtype, c_comm, &c_request);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
    return;
}

}


int
MPI_Alltoallv(const void *sendbuf, const int *sendcounts,
              const int *sdispls, MPI_Datatype sendtype, void *recvbuf,
              const int *recvcounts, const int *rdispls, MPI_Datatype recvtype,
              MPI_Comm comm)
{
    int ret,sum,i,sz;
    double t_elapsed;
    sum = 0;
    t_elapsed = MPI_Wtime();
    if ( prof_enabled == 1 ){
        ret = PMPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts,
                             rdispls, recvtype, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;

        MPI_Comm_size(comm, &sz);
        for ( i=0; i<sz; i++ ){
            if ( sendcounts[i] > 0 )
                sum+=sendcounts[i];
        }
        /* We won't need this reduce just sum all in the end */
        /* PMPI_Reduce(&sum, &sum_max, 1, MPI_INT, MPI_MAX, 0, comm); */
        profile_this(comm,sum,sendtype,Alltoallv,t_elapsed,1);
    }
    else{
        ret = PMPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf,
                             recvcounts, rdispls, recvtype, comm);
    }
    return ret;
}

extern "C" {
void
mpi_alltoallv_(void  *sendbuf, const int  *sendcnts, const int  *sdispls,
                  MPI_Fint  * sendtype, void  *recvbuf, const int  *recvcnts,
                  const int  *rdispls, MPI_Fint  * recvtype, MPI_Fint  * comm ,
                  MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);


   // All processes must use MPI_IN_PLACE
    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Alltoallv(sendbuf, sendcnts, sdispls, c_sendtype, recvbuf,
                            recvcnts, rdispls, c_recvtype, c_comm);

    *ierr = ret;
}
}

int MPI_Ialltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                   MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                   const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                   MPI_Request *request)
{
    int ret,sum,i,sz;
    double t_elapsed;
    sum = 0;
    t_elapsed = MPI_Wtime();
    if ( prof_enabled == 1 ){
        ret = PMPI_Ialltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf,
                              recvcounts, rdispls, recvtype, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;

        MPI_Comm_size(comm, &sz);
        for ( i=0; i<sz; i++ ){
            if ( sendcounts[i] > 0 )
                sum+=sendcounts[i];
        }
        /* We won't need this reduce just sum all in the end */
        /* PMPI_Reduce(&sum, &sum_max, 1, MPI_INT, MPI_MAX, 0, comm); */
        profile_this(comm,sum,sendtype,Ialltoallv,t_elapsed,1);
        requests_map[*request] = comm;
    }
    else{
        ret = PMPI_Ialltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf,
                              recvcounts, rdispls, recvtype, comm, request);
    }
    return ret;
}

extern "C" {
void
mpi_ialltoallv_(void *sendbuf, const int *sendcnts, const int *sdispls,
                            MPI_Fint *sendtype, void *recvbuf, const int *recvcnts,
                            const int *rdispls, MPI_Fint *recvtype, MPI_Fint *comm,
                            MPI_Fint *request, MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);


   // All processes must use MPI_IN_PLACE
    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Ialltoallv(sendbuf, sendcnts, sdispls, c_sendtype, recvbuf,
                         recvcnts, rdispls, c_recvtype, c_comm, &c_request);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);

}
}


int
MPI_Alltoallw(const void *sendbuf, const int *sendcounts, const int *sdispls,
              const MPI_Datatype *sendtypes, void *recvbuf, const int *recvcounts,
              const int *rdispls, const MPI_Datatype *recvtypes, MPI_Comm comm)
{
    int ret,sum,i,sz,type_sz;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        sum = 0;
        t_elapsed = MPI_Wtime();
        ret = PMPI_Alltoallw(sendbuf, sendcounts, sdispls, sendtypes, recvbuf,
                             recvcounts, rdispls, recvtypes, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;

        MPI_Comm_size(comm, &sz);
        for ( i=0; i<sz; i++ ){
            if ( recvcounts[i] > 0 ){
                PMPI_Type_size(recvtypes[i], &type_sz);
                sum+=(recvcounts[i]*type_sz);
            }
        }
        /* We won't need this reduce just sum all in the end */
        /* PMPI_Reduce(&sum, &sum_max, 1, MPI_INT, MPI_MAX, 0, comm); */
        profile_this(comm,sum,MPI_DATATYPE_NULL,Alltoallw,t_elapsed,1);
    }
    else{
        ret = PMPI_Alltoallw(sendbuf, sendcounts, sdispls, sendtypes, recvbuf,
                             recvcounts, rdispls, recvtypes, comm);
    }
    return ret;
}

extern "C" {
void mpi_alltoallw_(void *sendbuf, const int *sendcnts, const int *sdispls,
                        MPI_Fint * sendtypes, void *recvbuf, const int *recvcnts,
                        const int *rdispls, MPI_Fint *recvtypes, MPI_Fint *comm ,
                        MPI_Fint *ierr)
{
    int ret, comm_sz,i;
    MPI_Datatype *c_sendtypes = NULL;
    MPI_Datatype *c_recvtypes = NULL;
    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);
    PMPI_Comm_size(c_comm, &comm_sz);

    // All processes must use MPI_IN_PLACE
    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    c_recvtypes = (MPI_Datatype *)malloc(sizeof(MPI_Datatype) * comm_sz);
    if ( sendbuf != MPI_IN_PLACE){
        c_sendtypes = (MPI_Datatype *)malloc(sizeof(MPI_Datatype) * comm_sz);
        for (i=0; i<comm_sz; i++){
            c_sendtypes[i] = MPI_Type_f2c(sendtypes[i]);
            c_recvtypes[i] = MPI_Type_f2c(recvtypes[i]);
        }
    }
    else{
        for (i=0; i<comm_sz; i++){
            c_recvtypes[i] = MPI_Type_f2c(recvtypes[i]);
        }
    }

    ret = MPI_Alltoallw(sendbuf, sendcnts, sdispls, c_sendtypes, recvbuf,
                        recvcnts, rdispls, c_recvtypes, c_comm);
    *ierr = ret;

    if ( sendbuf != MPI_IN_PLACE)
        free(c_sendtypes);
    free(c_recvtypes);
    return;
}
}

int
MPI_Ialltoallw(const void *sendbuf, const int sendcounts[], const int sdispls[],
                   const MPI_Datatype sendtypes[], void *recvbuf,
                   const int recvcounts[], const int rdispls[],
                   const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request)
{
    int ret,sum,i,sz,type_sz;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        sum = 0;
        t_elapsed = MPI_Wtime();
        ret = PMPI_Ialltoallw(sendbuf, sendcounts, sdispls, sendtypes, recvbuf,
                              recvcounts, rdispls, recvtypes, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;

        MPI_Comm_size(comm, &sz);
        for ( i=0; i<sz; i++ ){
            if ( recvcounts[i] > 0 ){
                PMPI_Type_size(recvtypes[i], &type_sz);
                sum+=(recvcounts[i]*type_sz);
            }
        }
        /* We won't need this reduce just sum all in the end */
        /* PMPI_Reduce(&sum, &sum_max, 1, MPI_INT, MPI_MAX, 0, comm); */
        profile_this(comm,sum,MPI_DATATYPE_NULL,Ialltoallw,t_elapsed,1);
        requests_map[*request] = comm;
    }
    else{
        ret = PMPI_Ialltoallw(sendbuf, sendcounts, sdispls, sendtypes, recvbuf,
                              recvcounts, rdispls, recvtypes, comm, request);
    }
    return ret;
}

extern "C" {
void
mpi_ialltoallw_(const void *sendbuf, const int *sendcnts, const int *sdispls,
                        MPI_Fint *sendtypes, void *recvbuf, const int *recvcnts,
                        const int *rdispls, MPI_Fint *recvtypes, MPI_Fint *comm,
                        MPI_Fint *request, MPI_Fint *ierr)
{

    int ret, comm_sz,i;
    MPI_Datatype *c_sendtypes = NULL;
    MPI_Datatype *c_recvtypes = NULL;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_comm = MPI_Comm_f2c(*comm);
    PMPI_Comm_size(c_comm, &comm_sz);

    // All processes must use MPI_IN_PLACE
    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    c_recvtypes = (MPI_Datatype *)malloc(sizeof(MPI_Datatype) * comm_sz);
    if ( sendbuf != MPI_IN_PLACE){
        c_sendtypes = (MPI_Datatype *)malloc(sizeof(MPI_Datatype) * comm_sz);
        for (i=0; i<comm_sz; i++){
            c_sendtypes[i] = MPI_Type_f2c(sendtypes[i]);
            c_recvtypes[i] = MPI_Type_f2c(recvtypes[i]);
        }
    }
    else{
        for (i=0; i<comm_sz; i++){
            c_recvtypes[i] = MPI_Type_f2c(recvtypes[i]);
        }
    }

    ret = MPI_Ialltoallw(sendbuf, sendcnts, sdispls, c_sendtypes, recvbuf,
                        recvcnts, rdispls, c_recvtypes, c_comm, &c_request);

    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
    return;
}
}


int
MPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
               void *recvbuf, const int *recvcounts, const int *displs,
               MPI_Datatype recvtype, MPI_Comm comm)
{
    int ret,sum,i,sz;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts,
                              displs, recvtype, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;

        sum = 0;
        MPI_Comm_size(comm, &sz);
        for ( i=0; i<sz; i++ ){
            if( recvcounts[i] > 0 )
                sum+=recvcounts[i];
        }

        profile_this(comm,sum,recvtype,Allgatherv,t_elapsed,1);
    }
    else{
        ret = PMPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts,
                              displs, recvtype, comm);
    }
    return ret;
}


extern "C" {
void
mpi_allgatherv_(void  *sendbuf, int  * sendcount, MPI_Fint  * sendtype,
                   void  *recvbuf, int  *recvcounts, int  *displs,
                   MPI_Fint  * recvtype, MPI_Fint  * comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);


    // All processes must use MPI_IN_PLACE
    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Allgatherv(sendbuf, *sendcount, c_sendtype, recvbuf, recvcounts, displs, c_recvtype, c_comm);

    *ierr = ret;
    return;
}
}

int
MPI_Iallgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                void *recvbuf, const int *recvcounts, const int *displs,
                MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
{
    int ret,sum,i,sz;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Iallgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts,
                               displs, recvtype, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;

        sum = 0;
        MPI_Comm_size(comm, &sz);
        for ( i=0; i<sz; i++ ){
            if( recvcounts[i] > 0 )
                sum+=recvcounts[i];
        }

        profile_this(comm,sum,recvtype,Iallgatherv,t_elapsed,1);
        requests_map[*request] = comm;
    }
    else{
        ret = PMPI_Iallgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts,
                               displs, recvtype, comm, request);
    }
    return ret;

}

extern "C" {
void
mpi_iallgatherv_(void *sendbuf, int *sendcount, MPI_Fint *sendtype,
                    void  *recvbuf, int *recvcounts, int *displs,
                    MPI_Fint  *recvtype, MPI_Fint  *comm, MPI_Fint *request, MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);


    // All processes must use MPI_IN_PLACE
    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Iallgatherv(sendbuf, *sendcount, c_sendtype, recvbuf,
                          recvcounts, displs, c_recvtype, c_comm, &c_request);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
}
}


int
MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
           MPI_Op op, int root, MPI_Comm comm)
{
    int ret;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();

        ret = PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;

        profile_this(comm,count,datatype,Reduce,t_elapsed,0);
    }
    else{
        ret = PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
    }
    return ret;
}

extern "C" {
void
mpi_reduce_(void  *sendbuf, void  *recvbuf, int  * count,
                    MPI_Fint  * datatype, MPI_Fint  * op, int  * root,
                    MPI_Fint  * comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Op c_op;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_op = MPI_Op_f2c(*op);
    c_comm = MPI_Comm_f2c(*comm);


    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Reduce(sendbuf, recvbuf, *count, c_datatype, c_op, *root, c_comm);

    *ierr = ret;
    return;

}
}

int
MPI_Ireduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                MPI_Op op, int root, MPI_Comm comm, MPI_Request *request)
{
    int ret;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Ireduce(sendbuf, recvbuf, count, datatype, op, root, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;

        profile_this(comm,count,datatype,Ireduce,t_elapsed,0);
        requests_map[*request] = comm;
    }
    else{
        ret = PMPI_Ireduce(sendbuf, recvbuf, count, datatype, op, root, comm, request);
    }
    return ret;
}

extern "C" {
void
mpi_ireduce_(void  *sendbuf, void  *recvbuf, int  * count,
                     MPI_Fint  * datatype, MPI_Fint  * op, int  * root,
                     MPI_Fint  * comm , MPI_Fint  *request , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Op c_op;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_datatype = MPI_Type_f2c(*datatype);
    c_op = MPI_Op_f2c(*op);
    c_comm = MPI_Comm_f2c(*comm);


    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Ireduce(sendbuf, recvbuf, *count, c_datatype, c_op, *root, c_comm, &c_request);

    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
    return;
}
}

int
MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
           int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    int ret;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);

        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm,sendcount,sendtype,Gather,t_elapsed,0);
    }
    else{
        ret = PMPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
    }
    return ret;

}

extern "C" {
void
mpi_gather_(const void  *sendbuf, int  * sendcnt, MPI_Fint  * sendtype,
               void  *recvbuf, int  * recvcount, MPI_Fint  * recvtype,
               int  * root, MPI_Fint  * comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    // Only the root process can use MPI_IN_PLACE
    if ( sendbuf == mpisee_fortran_mpi_in_place)
        ret = MPI_Gather(MPI_IN_PLACE, *sendcnt, c_sendtype, recvbuf, *recvcount, c_recvtype, *root, c_comm);
    else
        ret = MPI_Gather(sendbuf, *sendcnt, c_sendtype, recvbuf, *recvcount, c_recvtype, *root, c_comm);

    *ierr = ret;
    return;

}
}

int
MPI_Igather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                void *recvbuf, int recvcount, MPI_Datatype recvtype,
                int root, MPI_Comm comm, MPI_Request *request)
{
    int ret;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Igather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, request);

        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm,sendcount,sendtype,Igather,t_elapsed,0);
        requests_map[*request] = comm;
    }
    else{
        ret = PMPI_Igather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, request);
    }
    return ret;

}

extern "C" {
void mpi_igather_(void  *sendbuf, int  * sendcnt, MPI_Fint  * sendtype,
                     void  *recvbuf, int  * recvcount, MPI_Fint  * recvtype,
                     int  * root, MPI_Fint  * comm , MPI_Fint  *request , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Igather(sendbuf, *sendcnt, c_sendtype, recvbuf, *recvcount, c_recvtype, *root, c_comm, &c_request);

    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
    return;
}
}

int
MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
            void *recvbuf, const int *recvcounts, const int *displs,
            MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    int ret,sum,rank,comm_size;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs,
                           recvtype, root, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;

        sum = 0;

        PMPI_Comm_rank(comm, &rank);
        if ( rank == root ){
            PMPI_Comm_size(comm, &comm_size);
            for (int i = 0; i < comm_size; i++){
                sum += recvcounts[i];
            }
        }
        else{
            sum = 0;
        }
        profile_this(comm, sum, recvtype, Gatherv, t_elapsed, 1);
    }
    else{
        ret = PMPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs,
                           recvtype, root, comm);
    }

    return ret;
}



extern "C" {
void
mpi_gatherv_(void  *sendbuf, int* sendcount, MPI_Fint  * sendtype,
                void  *recvbuf, const int* recvcounts, const int  *displs,
                MPI_Fint  * recvtype, int  * root, MPI_Fint  * comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Gatherv(sendbuf, *sendcount, c_sendtype, recvbuf, recvcounts,
                      displs, c_recvtype, *root, c_comm);

    *ierr = (MPI_Fint)ret;
    return;
}
}

int
MPI_Igatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
             void *recvbuf, const int *recvcounts, const int *displs,
             MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request)
{
    int ret,sum,rank,comm_size;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Igatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs,
                            recvtype, root, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;

        sum = 0;

        PMPI_Comm_rank(comm, &rank);
        if ( rank == root ){
            PMPI_Comm_size(comm, &comm_size);
            for (int i = 0; i < comm_size; i++){
                sum += recvcounts[i];
            }
        }
        else{
            sum = 0;
        }
        profile_this(comm, sum, recvtype, Igatherv, t_elapsed, 1);
        requests_map[*request] = comm;
    }
    else{
        ret = PMPI_Igatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs,
                            recvtype, root, comm, request);
    }
    return ret;
}

extern "C" {
void
mpi_igatherv_(void  *sendbuf, int  * sendcount, MPI_Fint  * sendtype,
                 void  *recvbuf, const int  *recvcounts, const int  *displs,
                 MPI_Fint  * recvtype, int  * root, MPI_Fint  * comm , MPI_Fint  *request , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Igatherv(sendbuf, *sendcount, c_sendtype, recvbuf, recvcounts, displs,
                           c_recvtype, *root, c_comm, &c_request);

    *ierr = (MPI_Fint)ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
    return;
}
}


int
MPI_Scatterv(const void *sendbuf, const int sendcounts[], const int *displs,
             MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, int root, MPI_Comm comm)
{

    int ret,rank,comm_size;
    uint64_t sum;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        sum = 0;
        t_elapsed = MPI_Wtime();
        ret = PMPI_Scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount,
                            recvtype, root, comm);

        t_elapsed = MPI_Wtime() - t_elapsed;
        PMPI_Comm_rank(comm, &rank);
        if ( rank == root ){
            PMPI_Comm_size(comm, &comm_size);
            for (int i = 0; i < comm_size; i++){
                sum += sendcounts[i];
            }
        }
        else{
            sum = 0;
        }
        profile_this(comm, sum, sendtype, Scatterv, t_elapsed, 1);
    }
    else{
        ret = PMPI_Scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount,
                            recvtype, root, comm);
    }
    return ret;
}

extern "C" {
void
mpi_scatterv_(void  *sendbuf, const int  *sendcounts, const int  *displs,
                 MPI_Fint  * sendtype, void  *recvbuf, int  * recvcount,
                 MPI_Fint  * recvtype, int  * root, MPI_Fint  * comm, MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Scatterv(sendbuf, sendcounts, displs, c_sendtype, recvbuf, *recvcount,
                 c_recvtype, *root, c_comm);

    *ierr = (MPI_Fint)ret;
    return;
}
}

int
MPI_Iscatterv(const void *sendbuf, const int sendcounts[], const int displs[],
                  MPI_Datatype sendtype, void *recvbuf, int recvcount,
                  MPI_Datatype recvtype, int root, MPI_Comm comm, MPI_Request *request)
{
    int ret,rank,comm_size;
    uint64_t sum;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        sum = 0;
        t_elapsed = MPI_Wtime();
        ret = PMPI_Iscatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount,
                             recvtype, root, comm, request);

        t_elapsed = MPI_Wtime() - t_elapsed;
        PMPI_Comm_rank(comm, &rank);
        if ( rank == root ){
            PMPI_Comm_size(comm, &comm_size);
            for (int i = 0; i < comm_size; i++){
                sum += sendcounts[i];
            }
        }
        else{
            sum = 0;
        }
        profile_this(comm, sum, sendtype, Iscatterv, t_elapsed, 1);
        requests_map[*request] = comm;
    }
    else{
        ret = PMPI_Iscatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount,
                             recvtype, root, comm, request);
    }
    return ret;
}

extern "C" {
void
mpi_iscatterv_(void  *sendbuf, const int  *sendcounts, const int  *displs,
                  MPI_Fint  * sendtype, void  *recvbuf, int  * recvcount,
                  MPI_Fint  * recvtype, int  * root, MPI_Fint  * comm, MPI_Fint  *request, MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Iscatterv(sendbuf, sendcounts, displs, c_sendtype, recvbuf, *recvcount,
                        c_recvtype, *root, c_comm, &c_request);

    *ierr = (MPI_Fint)ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
    return;
}
}

int
MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
            void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
            MPI_Comm comm)
{
    int ret,rank,sum;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                           recvtype, root, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;

        PMPI_Comm_rank(comm, &rank);
        if ( rank == root ){
            sum = sendcount;
        }
        else{
            sum = 0;
        }
        profile_this(comm, sum, sendtype, Scatter, t_elapsed, 0);
    }
    else{
        ret = PMPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                           recvtype, root, comm);
    }
    return ret;
}

extern "C" {
void
mpi_scatter_(void  *sendbuf, int  * sendcount, MPI_Fint  * sendtype,
                void  *recvbuf, int  * recvcount, MPI_Fint  * recvtype, int  * root,
                MPI_Fint  * comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Scatter(sendbuf, *sendcount, c_sendtype, recvbuf, *recvcount,
                      c_recvtype, *root, c_comm);

    *ierr = (MPI_Fint)ret;
    return;
}
}

int
MPI_Iscatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
             void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
             MPI_Comm comm, MPI_Request *request)
{
    int ret,rank,sum;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Iscatter(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                            recvtype, root, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;

        PMPI_Comm_rank(comm, &rank);
        if ( rank == root ){
            sum = sendcount;
        }
        else{
            sum = 0;
        }
        profile_this(comm, sum, sendtype, Iscatter, t_elapsed, 0);
        requests_map[*request] = comm;
    }
    else{
        ret = PMPI_Iscatter(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                            recvtype, root, comm, request);
    }
    return ret;
}

extern "C" {
void
mpi_iscatter_(void  *sendbuf, int  * sendcount, MPI_Fint  * sendtype,
                 void  *recvbuf, int  * recvcount, MPI_Fint  * recvtype, int  * root,
                 MPI_Fint  * comm , MPI_Fint  *request , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Iscatter(sendbuf, *sendcount, c_sendtype, recvbuf, *recvcount,
                        c_recvtype, *root, c_comm, &c_request);

    *ierr = (MPI_Fint)ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
    return;
}
}


int
MPI_Scan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
         MPI_Op op, MPI_Comm comm)
{
    int ret;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Scan(sendbuf, recvbuf, count, datatype, op, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;

        profile_this(comm,count,datatype,Scan,t_elapsed,0);
    }
    else{
        ret = PMPI_Scan(sendbuf, recvbuf, count, datatype, op, comm);
    }
    return ret;

}


extern "C" {
void
mpi_scan_(void  *sendbuf, void  *recvbuf, int  * count, MPI_Fint  * datatype,
             MPI_Fint  * op, MPI_Fint  * comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Op c_op;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_op = PMPI_Op_f2c(*op);
    c_comm = MPI_Comm_f2c(*comm);

    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Scan(sendbuf, recvbuf, *count, c_datatype, c_op, c_comm);


    *ierr = (MPI_Fint)ret;
    return;
}
}

int
MPI_Iscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
              MPI_Op op, MPI_Comm comm, MPI_Request *request)
{
    int ret;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Iscan(sendbuf, recvbuf, count, datatype, op, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;

        profile_this(comm,count,datatype,Iscan,t_elapsed,0);
        requests_map[*request] = comm;
    }
    else{
        ret = PMPI_Iscan(sendbuf, recvbuf, count, datatype, op, comm, request);
    }
    return ret;
}

extern "C" {
void
mpi_iscan_(void  *sendbuf, void  *recvbuf, int  * count, MPI_Fint  * datatype,
              MPI_Fint  * op, MPI_Fint  * comm , MPI_Fint  *request , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Op c_op;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_datatype = MPI_Type_f2c(*datatype);
    c_op = PMPI_Op_f2c(*op);
    c_comm = MPI_Comm_f2c(*comm);

    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Iscan(sendbuf, recvbuf, *count, c_datatype, c_op, c_comm, &c_request);

    *ierr = (MPI_Fint)ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
    return;
}
}

int
MPI_Exscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
           MPI_Op op, MPI_Comm comm)
{
    int ret;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Exscan(sendbuf, recvbuf, count, datatype, op, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm,count,datatype,Exscan,t_elapsed,0);
    }
    else{
        ret = PMPI_Exscan(sendbuf, recvbuf, count, datatype, op, comm);
    }
    return ret;
}

extern "C"{
void
mpi_exscan_(void  *sendbuf, void  *recvbuf, int  * count, MPI_Fint  * datatype,
               MPI_Fint  * op, MPI_Fint  * comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Op c_op;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_op = MPI_Op_f2c(*op);
    c_comm = MPI_Comm_f2c(*comm);

    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    if ( sendbuf == mpisee_fortran_mpi_in_place)
        ret = MPI_Exscan(MPI_IN_PLACE, recvbuf, *count, c_datatype, c_op, c_comm);
    else
        ret = MPI_Exscan(sendbuf, recvbuf, *count, c_datatype, c_op, c_comm);

    *ierr = ret;
    return;
}
}

int
MPI_Iexscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
            MPI_Op op, MPI_Comm comm, MPI_Request *request)
{
    int ret;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Iexscan(sendbuf, recvbuf, count, datatype, op, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm,count,datatype,Iexscan,t_elapsed,0);
        requests_map[*request] = comm;
    }
    else{
        ret = PMPI_Iexscan(sendbuf, recvbuf, count, datatype, op, comm, request);
    }
    return ret;
}

extern "C"{
void
mpi_iexscan_(void  *sendbuf, void  *recvbuf, int  * count, MPI_Fint  * datatype,
                MPI_Fint  * op, MPI_Fint  * comm, MPI_Fint  *request, MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Op c_op;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_datatype = MPI_Type_f2c(*datatype);
    c_op = MPI_Op_f2c(*op);
    c_comm = MPI_Comm_f2c(*comm);

    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Iexscan(sendbuf, recvbuf, *count, c_datatype, c_op, c_comm, &c_request);

    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
    return;
}
}
int
MPI_Barrier ( MPI_Comm comm )
{
    int ret;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Barrier(comm);
        t_elapsed = MPI_Wtime()-t_elapsed;

        profile_this(comm,0,MPI_DATATYPE_NULL,Barrier,t_elapsed,1);
    }
    else{
        ret = PMPI_Barrier(comm);
    }
    return ret;
}


extern "C" {
void
mpi_barrier_(MPI_Fint  * comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Barrier(c_comm);
    *ierr = (MPI_Fint)ret;
    return;
}
}

int
MPI_Ibarrier(MPI_Comm comm, MPI_Request *request)
{
    int ret;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Ibarrier(comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm,0,MPI_DATATYPE_NULL,Ibarrier,t_elapsed,1);
        requests_map[*request] = comm;
    }
    else{
        ret = PMPI_Ibarrier(comm, request);
    }
    return ret;
}

extern "C"{
void
mpi_ibarrier_(MPI_Fint *comm, MPI_Fint *request, MPI_Fint *ierr)
{
    int ret;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Ibarrier(c_comm, &c_request);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
    return;
}
}

int
MPI_Reduce_scatter(const void *sendbuf, void *recvbuf, const int *recvcounts,
                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
    int ret,rank,sum;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Reduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm);
        t_elapsed = MPI_Wtick() - t_elapsed;
        PMPI_Comm_rank(comm, &rank);
        sum = recvcounts[rank];
        profile_this(comm, sum, datatype, Reduce_scatter, t_elapsed, 1);
    }
    else{
        ret = PMPI_Reduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm);
    }
    return ret;

}

extern "C" {
void
mpi_reduce_scatter_(void  *sendbuf, void  *recvbuf, const int *recvcnts,
                            MPI_Fint  * datatype, MPI_Fint  * op, MPI_Fint  * comm,
                            MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Op c_op;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_op = MPI_Op_f2c(*op);
    c_comm = MPI_Comm_f2c(*comm);

    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Reduce_scatter(sendbuf, recvbuf, recvcnts, c_datatype, c_op, c_comm);

    *ierr = (MPI_Fint)ret;
    return;

}
}

int
MPI_Ireduce_scatter(const void *sendbuf, void *recvbuf, const int *recvcounts,
                    MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request)
{
    int ret,sum,rank;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Ireduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;
        PMPI_Comm_rank(comm, &rank);
        sum = recvcounts[rank];
        profile_this(comm, sum, datatype, Ireduce_scatter, t_elapsed, 1);
        requests_map[*request] = comm;
    }
    else{
        ret = PMPI_Ireduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm, request);
    }
    return ret;
}

extern "C" {
void
mpi_ireduce_scatter_(void  *sendbuf, void  *recvbuf, const int *recvcnts,
                            MPI_Fint  * datatype, MPI_Fint  * op, MPI_Fint  * comm,
                            MPI_Fint  *request , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Op c_op;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_datatype = MPI_Type_f2c(*datatype);
    c_op = MPI_Op_f2c(*op);
    c_comm = MPI_Comm_f2c(*comm);

    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Ireduce_scatter(sendbuf, recvbuf, recvcnts, c_datatype, c_op, c_comm, &c_request);

    *ierr = (MPI_Fint)ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
    return;

}
}

int
MPI_Reduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount,
                         MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
    int ret;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Reduce_scatter_block(sendbuf, recvbuf, recvcount, datatype, op, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm, recvcount, datatype, Reduce_scatter_block, t_elapsed, 0);
    }
    else{
        ret = PMPI_Reduce_scatter_block(sendbuf, recvbuf, recvcount, datatype, op, comm);
    }
    return ret;

}

extern "C" {
void
mpi_reduce_scatter_block_(void *sendbuf, void *recvbuf, int *recvcount,
                            MPI_Fint *datatype, MPI_Fint *op, MPI_Fint *comm,
                            MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Op c_op;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_op = MPI_Op_f2c(*op);
    c_comm = MPI_Comm_f2c(*comm);

    if ( sendbuf == mpisee_fortran_mpi_in_place)
        sendbuf = MPI_IN_PLACE;
    if ( sendbuf == mpisee_fortran_mpi_bottom )
        sendbuf = MPI_BOTTOM;
    if ( recvbuf == mpisee_fortran_mpi_bottom )
        recvbuf = MPI_BOTTOM;

    ret = MPI_Reduce_scatter_block(sendbuf, recvbuf, *recvcount, c_datatype, c_op, c_comm);

    *ierr = (MPI_Fint)ret;
    return;

}
}

// Missing Ireduce_scatter_block
