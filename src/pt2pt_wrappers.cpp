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
#include "commprof.h"
#include <mpi.h>

int
MPI_Isend(const void *buf, int count, MPI_Datatype datatype,int dest, int tag,
          MPI_Comm comm, MPI_Request *request)
{
    int ret;
    double t_elapsed;
    if ( prof_enabled  == 1){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Isend(buf, count, datatype, dest, tag, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm, count, datatype, Isend, t_elapsed, 0);
        requests_map[*request] = comm;
    }
    else{
        ret = PMPI_Isend(buf, count, datatype, dest, tag, comm, request);
    }
    return ret;
}

extern "C" {
void
mpi_isend_(const void  *buf, int  * count, MPI_Fint  * datatype,
                          int  * dest, int  * tag, MPI_Fint  * comm,
                          MPI_Fint  *request , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    MPI_Request c_request;
    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Isend(buf, *count, c_datatype, *dest, *tag, c_comm, &c_request);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
    return;
}
}


extern "C" {
int
MPI_Send(const void *buf, int count, MPI_Datatype datatype,
         int dest,int tag, MPI_Comm comm)
{
    int ret;
    double t_elapsed = 0.0;
    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Send(buf, count, datatype, dest, tag, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm, count, datatype, Send, t_elapsed, 0);
    }
    else{
        ret = PMPI_Send(buf, count, datatype, dest, tag, comm);
    }
    return ret;
}
}

extern "C" {
void
mpi_send_(const void  *buf, int  * count, MPI_Fint  * datatype,
             int  * dest, int  * tag, MPI_Fint  * comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Send(buf, *count, c_datatype, *dest, *tag, c_comm);
    *ierr = ret;
    return;
}
}

extern "C"{
int
MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
          MPI_Comm comm, MPI_Request *request)
{

    int ret;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Irecv(buf, count, datatype, source, tag, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm, count,datatype,Irecv,t_elapsed,0);
        requests_map[*request] = comm;
    }
    else{
        ret = PMPI_Irecv(buf, count, datatype, source, tag, comm, request);
    }
    return ret;
}
}

extern "C" {
void
mpi_irecv_(void  *buf, int  * count, MPI_Fint  * datatype, int  * source,
              int  * tag, MPI_Fint  * comm, MPI_Fint  *request , MPI_Fint *ierr)
{

    int ret;
    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Irecv(buf, *count, c_datatype, *source, *tag, c_comm, &c_request);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
    return;
}
}

int
MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
             MPI_Comm comm, MPI_Status *status)
{
    int ret;
    double t_elapsed;
    if ( prof_enabled == 1){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Recv(buf,count,datatype,source,tag,comm,status);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm,count,datatype,Recv,t_elapsed,0);
    }
    else{
        ret = PMPI_Recv(buf,count,datatype,source,tag,comm,status);
    }
    return ret;
}

extern "C" {
void
mpi_recv_(void* buf, int* count,MPI_Fint* datatype, int* source, int* tag,
             MPI_Fint  * comm, MPI_Status  *status , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Recv(buf,*count,c_datatype,*source,*tag,c_comm,status);
    *ierr = ret;
    return;
}
}

int
MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
             int dest, int sendtag, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, int source, int recvtag,
             MPI_Comm comm, MPI_Status *status)
{

    int ret;
    double t_elapsed;
    int64_t sum;

    if ( prof_enabled == 1){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf,
                            recvcount, recvtype, source, recvtag, comm, status);
        t_elapsed = MPI_Wtime() - t_elapsed;
        if ( dest == MPI_PROC_NULL )
            sum = 0;
        else
            sum = sendcount;
        // sum = sum | 0x1;
        // sum = sum>>1;
        profile_this(comm,sum,sendtype,Sendrecv,t_elapsed,0);
    }
    else{
        ret = PMPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf,
                        recvcount, recvtype, source, recvtag, comm, status);
    }

    return ret;
}


extern "C" {
void
mpi_sendrecv_(const void  *sendbuf, int  * sendcount,MPI_Fint  * sendtype,
                 int  * dest, int  * sendtag, void  *recvbuf, int  * recvcount,
                 MPI_Fint  * recvtype, int  * source, int  * recvtag,
                 MPI_Fint  * comm, MPI_Status  *status , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;
    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Sendrecv(sendbuf, *sendcount, c_sendtype, *dest, *sendtag, recvbuf,
                        *recvcount, c_recvtype, *source, *recvtag, c_comm, status);
    *ierr = ret;
    return;
}
}

// int
// MPI_Isendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
//                   int dest, int sendtag, void *recvbuf, int recvcount,
//                   MPI_Datatype recvtype, int source, int recvtag,
//                   MPI_Comm comm, MPI_Request *request)
// {
//     int ret;
//     double t_elapsed;
//     if ( prof_enabled == 1){
//         t_elapsed =  MPI_Wtime();
//         ret = PMPI_Isendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf,
//                              recvcount, recvtype, source, recvtag, comm, request);
//         t_elapsed = MPI_Wtime() - t_elapsed;
//         profile_this(comm,sendcount,sendtype,Isendrecv,t_elapsed,source);
//         requests_map[*request] = comm;
//     }
//     else{
//         ret = PMPI_Isendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf,
//                              recvcount, recvtype, source, recvtag, comm, request);
//     }
//     return ret;
// }

// extern "C" {
// void mpi_isendrecv_(const void  *sendbuf, int  * sendcount,MPI_Fint  * sendtype,
//                        int  * dest, int  * sendtag, void  *recvbuf, int  * recvcount,
//                        MPI_Fint  * recvtype, int  * source, int  * recvtag,
//                        MPI_Fint  * comm, MPI_Fint  *request , MPI_Fint *ierr)
// {
//     int ret;
//     MPI_Datatype c_sendtype;
//     MPI_Datatype c_recvtype;
//     MPI_Comm c_comm;
//     MPI_Request c_request;
//     c_sendtype = MPI_Type_f2c(*sendtype);
//     c_recvtype = MPI_Type_f2c(*recvtype);
//     c_comm = MPI_Comm_f2c(*comm);
//     ret = MPI_Isendrecv(sendbuf, *sendcount, c_sendtype, *dest, *sendtag, recvbuf,
//                         *recvcount, c_recvtype, *source, *recvtag, c_comm, &c_request);
//     *ierr = ret;
//     if ( ret == MPI_SUCCESS )
//         *request = MPI_Request_c2f(c_request);
//     return;
// }
// }


int
MPI_Ssend(const void *buf, int count, MPI_Datatype datatype, int dest,
              int tag, MPI_Comm comm)
{
    int ret;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Ssend(buf, count, datatype, dest, tag, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm,count,datatype,Ssend,t_elapsed,0);
    }
    else{
        ret = PMPI_Ssend(buf, count, datatype, dest, tag, comm);
    }
    return ret;
}

extern "C" {
void
mpi_ssend_(const void *buf, int *count, MPI_Fint  *datatype,
              int  *dest, int  *tag, MPI_Fint  *comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Ssend(buf, *count, c_datatype, *dest, *tag, c_comm);
    *ierr = ret;
}
}

int
MPI_Issend(const void *buf, int count, MPI_Datatype datatype, int dest,
               int tag, MPI_Comm comm, MPI_Request *request)
{
    int ret;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Issend(buf, count, datatype, dest, tag, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm, count, datatype, Issend, t_elapsed, 0);
        requests_map[*request] = comm;
    }
    else{
        ret = PMPI_Issend(buf, count, datatype, dest, tag, comm, request);
    }
    return ret;

}

extern "C"{
void
mpi_issend_(const void *buf, int *count, MPI_Fint *datatype,
               int *dest, int *tag, MPI_Fint *comm,
               MPI_Fint *request , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    MPI_Request c_request;
    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Issend(buf, *count, c_datatype, *dest, *tag, c_comm, &c_request);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
    return;
}
}

int
MPI_Bsend(const void *buf, int count, MPI_Datatype datatype, int dest,
               int tag, MPI_Comm comm)
{
    int ret;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Bsend(buf, count, datatype, dest, tag, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm, count, datatype, Bsend, t_elapsed, 0);
    }
    else{
        ret = PMPI_Bsend(buf, count, datatype, dest, tag, comm);
    }
    return ret;
}

extern "C" {
void
mpi_bsend_(const void *buf, int *count, MPI_Fint *datatype,
              int *dest, int *tag, MPI_Fint *comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Bsend(buf, *count, c_datatype, *dest, *tag, c_comm);
    *ierr = ret;
}
}

int
MPI_Ibsend(const void *buf, int count, MPI_Datatype datatype, int dest,
               int tag, MPI_Comm comm, MPI_Request *request)
{
    int ret;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Ibsend(buf, count, datatype, dest, tag, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm, count, datatype, Ibsend, t_elapsed, 0);
        requests_map[*request] = comm;
    }
    else{
        ret = PMPI_Ibsend(buf, count, datatype, dest, tag, comm, request);
    }
    return ret;
}

extern "C" {
void
mpi_ibsend_(const void *buf, int *count, MPI_Fint *datatype,
               int *dest, int *tag, MPI_Fint *comm,
               MPI_Fint *request , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    MPI_Request c_request;
    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Ibsend(buf, *count, c_datatype, *dest, *tag, c_comm, &c_request);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
    return;
}
}

int
MPI_Rsend(const void *buf, int count, MPI_Datatype datatype, int dest,
               int tag, MPI_Comm comm)
{
    int ret;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Rsend(buf, count, datatype, dest, tag, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm, count, datatype, Rsend, t_elapsed, 0);
    }
    else{
        ret = PMPI_Rsend(buf, count, datatype, dest, tag, comm);
    }
    return ret;
}

extern "C" {
void
mpi_rsend_(const void *buf, int *count, MPI_Fint *datatype,
               int *dest, int *tag, MPI_Fint *comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Rsend(buf, *count, c_datatype, *dest, *tag, c_comm);
    *ierr = ret;
}
}

int
MPI_Irsend(const void *buf, int count, MPI_Datatype datatype, int dest,
               int tag, MPI_Comm comm, MPI_Request *request)
{
    int ret;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Irsend(buf, count, datatype, dest, tag, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm, count, datatype, Irsend, t_elapsed, 0);
        requests_map[*request] = comm;
    }
    else{
        ret = PMPI_Irsend(buf, count, datatype, dest, tag, comm, request);
    }
    return ret;
}

extern "C" {
void
mpi_irsend_(const void *buf, int *count, MPI_Fint *datatype,
               int *dest, int *tag, MPI_Fint *comm,
               MPI_Fint *request , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    MPI_Request c_request;
    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Irsend(buf, *count, c_datatype, *dest, *tag, c_comm, &c_request);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
    return;
}
}
