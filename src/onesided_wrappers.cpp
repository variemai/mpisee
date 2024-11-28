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
#include "include/utils.h"
#include <mpi.h>


int
MPI_Win_create(void *base, MPI_Aint size, int disp_unit, MPI_Info info,
               MPI_Comm comm, MPI_Win *win)
{
    int ret;
    ret = PMPI_Win_create(base, size, disp_unit, info, comm, win);
    if (ret == MPI_SUCCESS)
    {
 #ifdef MPICH_NAME
        comm_map[*win] = comm;
 #endif
 #ifdef OMPI_MAJOR_VERSION
        MPI_Win_set_attr(*win, win_namekey(), comm);
 #endif


    }
    return ret;
}

extern "C"{
void mpi_win_create_(void *base, MPI_Aint *size, int *disp_unit, MPI_Fint *info,
                     MPI_Fint *comm, MPI_Fint *win, MPI_Fint *ierr)
{
    int ret;
    MPI_Info c_info = MPI_Info_f2c(*info);
    MPI_Comm c_comm = MPI_Comm_f2c(*comm);
    MPI_Win c_win;

    ret = MPI_Win_create(base, *size, *disp_unit, c_info, c_comm, &c_win);
    *win = MPI_Win_c2f(c_win);
    *ierr = ret;
}
}


int
MPI_Put(const void *origin_addr, int origin_count,
        MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
        int target_count, MPI_Datatype target_datatype, MPI_Win win)
{
    int ret;
    double t_elapsed;
    MPI_Comm comm = MPI_COMM_NULL;
    if (prof_enabled == 1){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Put(origin_addr, origin_count, origin_datatype, target_rank,
                   target_disp, target_count, target_datatype, win);
        t_elapsed = MPI_Wtime() - t_elapsed;
        #ifdef MPICH_NAME
        comm = comm_map[win];
        if (comm == MPI_COMM_NULL){
            mcpt_abort("Window communicator not found\n");
        }
        #endif
        #ifdef OMPI_MAJOR_VERSION
        int flag;
        MPI_Win_get_attr(win, win_namekey(), &comm, &flag);
        if (flag == 0){
            mcpt_abort("Window communicator not found\n");
        }
        #endif
        profile_this(comm, origin_count, origin_datatype, Put ,t_elapsed, 0);
    }
    else{
        ret = PMPI_Put(origin_addr, origin_count, origin_datatype, target_rank,
                   target_disp, target_count, target_datatype, win);
    }
    return ret;

}

extern "C" {
 void mpi_put_(const void *origin_addr, int *origin_count,
        MPI_Fint *origin_datatype, int *target_rank, int *target_disp,
        int *target_count, MPI_Fint *target_datatype, MPI_Fint *win, MPI_Fint *ierr)
{
    int ret;
    MPI_Win c_win = MPI_Win_f2c(*win);
    MPI_Datatype c_origin_datatype = MPI_Type_f2c(*origin_datatype);
    MPI_Datatype c_target_datatype = MPI_Type_f2c(*target_datatype);

    ret = MPI_Put(origin_addr, *origin_count, c_origin_datatype, *target_rank,
                  *target_disp, *target_count, c_target_datatype, c_win);
    *ierr = ret;
}
}

int
MPI_Rput(const void *origin_addr, int origin_count,
        MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
        int target_count, MPI_Datatype target_datatype, MPI_Win win, MPI_Request *request)
{
    int ret;
    double t_elapsed;
    MPI_Comm comm = MPI_COMM_NULL;
    if (prof_enabled == 1){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Rput(origin_addr, origin_count, origin_datatype, target_rank,
                   target_disp, target_count, target_datatype, win, request);
        t_elapsed = MPI_Wtime() - t_elapsed;
        #ifdef MPICH_NAME
        comm = comm_map[win];
        if (comm == MPI_COMM_NULL){
            mcpt_abort("Window communicator not found\n");
        }
        #endif
        #ifdef OMPI_MAJOR_VERSION
        int flag;
        MPI_Win_get_attr(win, win_namekey(), &comm, &flag);
        if (flag == 0){
            mcpt_abort("Window communicator not found\n");
        }
        #endif
        profile_this(comm, origin_count, origin_datatype, Rput ,t_elapsed, 0);
    }
    else{
        ret = PMPI_Rput(origin_addr, origin_count, origin_datatype, target_rank,
                   target_disp, target_count, target_datatype, win, request);
    }
    return ret;
}

extern "C" {
void mpi_rput_(const void *origin_addr, int *origin_count,
        MPI_Fint *origin_datatype, int *target_rank, int *target_disp,
        int *target_count, MPI_Fint *target_datatype, MPI_Fint *win, MPI_Fint *request, MPI_Fint *ierr)
{
    int ret;
    MPI_Win c_win = MPI_Win_f2c(*win);
    MPI_Datatype c_origin_datatype = MPI_Type_f2c(*origin_datatype);
    MPI_Datatype c_target_datatype = MPI_Type_f2c(*target_datatype);
    MPI_Request c_request;

    ret = MPI_Rput(origin_addr, *origin_count, c_origin_datatype, *target_rank,
                  *target_disp, *target_count, c_target_datatype, c_win, &c_request);
    *request = MPI_Request_c2f(c_request);
    *ierr = ret;
}
}

int
MPI_Accumulate(const void *origin_addr, int origin_count,
               MPI_Datatype origin_datatype, int target_rank,
               MPI_Aint target_disp, int target_count,
               MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
{
    int ret;
    double t_elapsed;
    MPI_Comm comm = MPI_COMM_NULL;
    if (prof_enabled == 1){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Accumulate(origin_addr, origin_count, origin_datatype, target_rank,
                   target_disp, target_count, target_datatype, op, win);
        t_elapsed = MPI_Wtime() - t_elapsed;
        #ifdef MPICH_NAME
        comm = comm_map[win];
        if (comm == MPI_COMM_NULL){
            mcpt_abort("Window communicator not found\n");
        }
        #endif
        #ifdef OMPI_MAJOR_VERSION
        int flag;
        MPI_Win_get_attr(win, win_namekey(), &comm, &flag);
        if (flag == 0){
            mcpt_abort("Window communicator not found\n");
        }
        #endif
        profile_this(comm, origin_count, origin_datatype, Accumulate ,t_elapsed, 0);
    }
    else{
        ret = PMPI_Accumulate(origin_addr, origin_count, origin_datatype, target_rank,
                   target_disp, target_count, target_datatype, op, win);
    }
    return ret;
}

extern "C" {
void mpi_accumulate_(const void *origin_addr, int *origin_count,
               MPI_Fint *origin_datatype, int *target_rank,
               int *target_disp, int *target_count,
               MPI_Fint *target_datatype, MPI_Fint *op, MPI_Fint *win, MPI_Fint *ierr)
{
    int ret;
    MPI_Win c_win = MPI_Win_f2c(*win);
    MPI_Datatype c_origin_datatype = MPI_Type_f2c(*origin_datatype);
    MPI_Datatype c_target_datatype = MPI_Type_f2c(*target_datatype);
    MPI_Op c_op = MPI_Op_f2c(*op);

    ret = MPI_Accumulate(origin_addr, *origin_count, c_origin_datatype, *target_rank,
                  *target_disp, *target_count, c_target_datatype, c_op, c_win);
    *ierr = ret;
}
}


int
MPI_Get(void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
        int target_rank, MPI_Aint target_disp, int target_count,
        MPI_Datatype target_datatype, MPI_Win win)
{
    int ret;
    double t_elapsed;
    MPI_Comm comm = MPI_COMM_NULL;
    if (prof_enabled == 1){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Get(origin_addr, origin_count, origin_datatype, target_rank,
                   target_disp, target_count, target_datatype, win);
        t_elapsed = MPI_Wtime() - t_elapsed;
        #ifdef MPICH_NAME
        comm = comm_map[win];
        if (comm == MPI_COMM_NULL){
            mcpt_abort("Window communicator not found\n");
        }
        #endif
        #ifdef OMPI_MAJOR_VERSION
        int flag;
        MPI_Win_get_attr(win, win_namekey(), &comm, &flag);
        if (flag == 0){
            mcpt_abort("Window communicator not found\n");
        }
        #endif
        profile_this(comm, target_count, target_datatype, Get ,t_elapsed, 0);
    }
    else{
        ret = PMPI_Get(origin_addr, origin_count, origin_datatype, target_rank,
                   target_disp, target_count, target_datatype, win);
    }
    return ret;
}

extern "C" {
void mpi_get_(void *origin_addr, int *origin_count, MPI_Fint *origin_datatype,
        int *target_rank, int *target_disp, int *target_count,
        MPI_Fint *target_datatype, MPI_Fint *win, MPI_Fint *ierr)
{
    int ret;
    MPI_Win c_win = MPI_Win_f2c(*win);
    MPI_Datatype c_origin_datatype = MPI_Type_f2c(*origin_datatype);
    MPI_Datatype c_target_datatype = MPI_Type_f2c(*target_datatype);

    ret = MPI_Get(origin_addr, *origin_count, c_origin_datatype, *target_rank,
                  *target_disp, *target_count, c_target_datatype, c_win);
    *ierr = ret;
}
}

int
MPI_Rget(void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
        int target_rank, MPI_Aint target_disp, int target_count,
        MPI_Datatype target_datatype, MPI_Win win, MPI_Request *request)
{
    int ret;
    double t_elapsed;
    MPI_Comm comm = MPI_COMM_NULL;
    if (prof_enabled == 1){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Rget(origin_addr, origin_count, origin_datatype, target_rank,
                   target_disp, target_count, target_datatype, win, request);
        t_elapsed = MPI_Wtime() - t_elapsed;
        #ifdef MPICH_NAME
        comm = comm_map[win];
        if (comm == MPI_COMM_NULL){
            mcpt_abort("Window communicator not found\n");
        }
        #endif
        #ifdef OMPI_MAJOR_VERSION
        int flag;
        MPI_Win_get_attr(win, win_namekey(), &comm, &flag);
        if (flag == 0){
            mcpt_abort("Window communicator not found\n");
        }
        #endif
        profile_this(comm, target_count, target_datatype, Rget ,t_elapsed, 0);
    }
    else{
        ret = PMPI_Rget(origin_addr, origin_count, origin_datatype, target_rank,
                   target_disp, target_count, target_datatype, win, request);
    }
    return ret;
}

extern "C" {
void mpi_rget_(void *origin_addr, int *origin_count, MPI_Fint *origin_datatype,
        int *target_rank, int *target_disp, int *target_count,
        MPI_Fint *target_datatype, MPI_Fint *win, MPI_Fint *request, MPI_Fint *ierr)
{
    int ret;
    MPI_Win c_win = MPI_Win_f2c(*win);
    MPI_Datatype c_origin_datatype = MPI_Type_f2c(*origin_datatype);
    MPI_Datatype c_target_datatype = MPI_Type_f2c(*target_datatype);
    MPI_Request c_request;

    ret = MPI_Rget(origin_addr, *origin_count, c_origin_datatype, *target_rank,
                  *target_disp, *target_count, c_target_datatype, c_win, &c_request);
    *request = MPI_Request_c2f(c_request);
    *ierr = ret;
}
}

int
MPI_Win_fence(int assert, MPI_Win win)
{
    int ret;
    double t_elapsed;
    MPI_Comm comm = MPI_COMM_NULL;
    if (prof_enabled == 1){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Win_fence(assert, win);
        t_elapsed = MPI_Wtime() - t_elapsed;
        #ifdef MPICH_NAME
        comm = comm_map[win];
        if (comm == MPI_COMM_NULL){
            mcpt_abort("Window communicator not found\n");
        }
        #endif
        #ifdef OMPI_MAJOR_VERSION
        int flag;
        MPI_Win_get_attr(win, win_namekey(), &comm, &flag);
        if (flag == 0){
            mcpt_abort("Window communicator not found\n");
        }
        #endif
        profile_this(comm, 0, MPI_DATATYPE_NULL, Fence ,t_elapsed, 1);
    }
    else{
        ret = PMPI_Win_fence(assert, win);
    }
    return ret;
}

extern "C" {
void mpi_win_fence_(int *assert, MPI_Fint *win, MPI_Fint *ierr)
{
    int ret;
    MPI_Win c_win = MPI_Win_f2c(*win);

    ret = MPI_Win_fence(*assert, c_win);
    *ierr = ret;
}
}

int
MPI_Win_start(MPI_Group group, int assert, MPI_Win win)
{
    int ret;
    double t_elapsed;
    MPI_Comm comm = MPI_COMM_NULL;
    if (prof_enabled == 1){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Win_start(group, assert, win);
        t_elapsed = MPI_Wtime() - t_elapsed;
        #ifdef MPICH_NAME
        comm = comm_map[win];
        if (comm == MPI_COMM_NULL){
            mcpt_abort("Window communicator not found\n");
        }
        #endif
        #ifdef OMPI_MAJOR_VERSION
        int flag;
        MPI_Win_get_attr(win, win_namekey(), &comm, &flag);
        if (flag == 0){
            mcpt_abort("Window communicator not found\n");
        }
        #endif
        profile_this(comm, 0, MPI_DATATYPE_NULL, Win_start ,t_elapsed, 1);
    }
    else{
        ret = PMPI_Win_start(group, assert, win);
    }
    return ret;
}

extern "C" {
void mpi_win_start_(MPI_Fint *group, int *assert, MPI_Fint *win, MPI_Fint *ierr)
{
    int ret;
    MPI_Win c_win = MPI_Win_f2c(*win);
    MPI_Group c_group = MPI_Group_f2c(*group);

    ret = MPI_Win_start(c_group, *assert, c_win);
    *ierr = ret;
}
}

int
MPI_Win_complete(MPI_Win win)
{
    int ret;
    double t_elapsed;
    MPI_Comm comm = MPI_COMM_NULL;
    if (prof_enabled == 1){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Win_complete(win);
        t_elapsed = MPI_Wtime() - t_elapsed;
        #ifdef MPICH_NAME
        comm = comm_map[win];
        if (comm == MPI_COMM_NULL){
            mcpt_abort("Window communicator not found\n");
        }
        #endif
        #ifdef OMPI_MAJOR_VERSION
        int flag;
        MPI_Win_get_attr(win, win_namekey(), &comm, &flag);
        if (flag == 0){
            mcpt_abort("Window communicator not found\n");
        }
        #endif
        profile_this(comm, 0, MPI_DATATYPE_NULL, Win_complete ,t_elapsed, 1);
    }
    else{
        ret = PMPI_Win_complete(win);
    }
    return ret;
}

extern "C" {
void mpi_win_complete_(MPI_Fint *win, MPI_Fint *ierr)
{
    int ret;
    MPI_Win c_win = MPI_Win_f2c(*win);

    ret = MPI_Win_complete(c_win);
    *ierr = ret;
}
}

int
MPI_Win_post(MPI_Group group, int assert, MPI_Win win)
{
    int ret;
    double t_elapsed;
    MPI_Comm comm = MPI_COMM_NULL;
    if (prof_enabled == 1){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Win_post(group, assert, win);
        t_elapsed = MPI_Wtime() - t_elapsed;
        #ifdef MPICH_NAME
        comm = comm_map[win];
        if (comm == MPI_COMM_NULL){
            mcpt_abort("Window communicator not found\n");
        }
        #endif
        #ifdef OMPI_MAJOR_VERSION
        int flag;
        MPI_Win_get_attr(win, win_namekey(), &comm, &flag);
        if (flag == 0){
            mcpt_abort("Window communicator not found\n");
        }
        #endif
        profile_this(comm, 0, MPI_DATATYPE_NULL, Win_post ,t_elapsed, 1);
    }
    else{
        ret = PMPI_Win_post(group, assert, win);
    }
    return ret;
}

extern "C" {
void mpi_win_post_(MPI_Fint *group, int *assert, MPI_Fint *win, MPI_Fint *ierr)
{
    int ret;
    MPI_Win c_win = MPI_Win_f2c(*win);
    MPI_Group c_group = MPI_Group_f2c(*group);

    ret = MPI_Win_post(c_group, *assert, c_win);
    *ierr = ret;
}
}

int
MPI_Win_wait(MPI_Win win)
{
    int ret;
    double t_elapsed;
    MPI_Comm comm = MPI_COMM_NULL;
    if (prof_enabled == 1){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Win_wait(win);
        t_elapsed = MPI_Wtime() - t_elapsed;
        #ifdef MPICH_NAME
        comm = comm_map[win];
        if (comm == MPI_COMM_NULL){
            mcpt_abort("Window communicator not found\n");
        }
        #endif
        #ifdef OMPI_MAJOR_VERSION
        int flag;
        MPI_Win_get_attr(win, win_namekey(), &comm, &flag);
        if (flag == 0){
            mcpt_abort("Window communicator not found\n");
        }
        #endif
        profile_this(comm, 0, MPI_DATATYPE_NULL, Win_wait ,t_elapsed, 1);
    }
    else{
        ret = PMPI_Win_wait(win);
    }
    return ret;
}

extern "C" {
void mpi_win_wait_(MPI_Fint *win, MPI_Fint *ierr)
{
    int ret;
    MPI_Win c_win = MPI_Win_f2c(*win);

    ret = MPI_Win_wait(c_win);
    *ierr = ret;
}
}

int
MPI_Win_test(MPI_Win win, int *flag)
{
    int ret, flag1;
    double t_elapsed;
    MPI_Comm comm;
    if (prof_enabled == 1){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Win_test(win, flag);
        t_elapsed = MPI_Wtime() - t_elapsed;
        MPI_Win_get_attr(win, win_namekey(), &comm, &flag1);
        if (flag1 == 0){
            mcpt_abort("Window communicator not found\n");
        }
        else{
            profile_this(comm, 0, MPI_DATATYPE_NULL, Win_test ,t_elapsed, 1);
        }
    }
    else{
        ret = PMPI_Win_test(win, flag);
    }
    return ret;
}

extern "C" {
void mpi_win_test_(MPI_Fint *win, int *flag, MPI_Fint *ierr)
{
    int ret;
    MPI_Win c_win = MPI_Win_f2c(*win);

    ret = MPI_Win_test(c_win, flag);
    *ierr = ret;
}
}
