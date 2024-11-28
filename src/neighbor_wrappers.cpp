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


#include <mpi.h>
#include "include/commprof.h"

int MPI_Neighbor_allgather(const void *sendbuf, int sendcount,
    MPI_Datatype sendtype, void *recvbuf, int recvcount,
    MPI_Datatype recvtype, MPI_Comm comm) {
    int ret;
    double t_elapsed;
    if ( prof_enabled == 1 ) {
        t_elapsed = MPI_Wtime();
        ret = PMPI_Neighbor_allgather(sendbuf, sendcount, sendtype, recvbuf,
                                      recvcount, recvtype, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm, sendcount, sendtype, Neighbor_allgather, t_elapsed,
            0);
    } else {
        ret = PMPI_Neighbor_allgather(sendbuf, sendcount, sendtype, recvbuf,
                                      recvcount, recvtype, comm);
    }
    return ret;
}

extern "C" {
void mpi_neighbor_allgather_(const void  *sendbuf, int  *sendcount,
    MPI_Fint  *sendtype, void  *recvbuf, int  *recvcount, MPI_Fint  *recvtype,
    MPI_Fint  *comm, MPI_Fint *ierr) {
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    // The IN_PLACE is not meaningful here

    ret = MPI_Neighbor_allgather(sendbuf, *sendcount, c_sendtype, recvbuf,
                                 *recvcount, c_recvtype, c_comm);
    *ierr = ret;
    return;
}
}

int MPI_Neighbor_allgatherv(const void *sendbuf, int sendcount,
    MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
    const int displs[], MPI_Datatype recvtype, MPI_Comm comm) {
    int ret;
    double t_elapsed;
    if ( prof_enabled == 1 ) {
        t_elapsed = MPI_Wtime();
        ret = PMPI_Neighbor_allgatherv(sendbuf, sendcount, sendtype, recvbuf,
            recvcounts, displs, recvtype, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm, sendcount, sendtype, Neighbor_allgatherv,
            t_elapsed, 1);
    } else {
        ret = PMPI_Neighbor_allgatherv(sendbuf, sendcount, sendtype, recvbuf,
            recvcounts, displs, recvtype, comm);
    }
    return ret;
}

extern "C" {
void mpi_neighbor_allgatherv_(const void  *sendbuf, int  *sendcount,
    MPI_Fint *sendtype, void  *recvbuf, const int  *recvcounts,
    const int *displs, MPI_Fint  *recvtype, MPI_Fint  *comm, MPI_Fint *ierr) {
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    // MPI_IN_PLACE is not meaningful here

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    ret = MPI_Neighbor_allgatherv(sendbuf, *sendcount, c_sendtype, recvbuf,
        recvcounts, displs, c_recvtype, c_comm);
    *ierr = ret;
}
}

int MPI_Neighbor_alltoall(const void *sendbuf, int sendcount,
    MPI_Datatype sendtype, void *recvbuf, int recvcount,
    MPI_Datatype recvtype, MPI_Comm comm) {
    int ret;
    double t_elapsed;
    if ( prof_enabled == 1 ) {
        t_elapsed = MPI_Wtime();
        ret = PMPI_Neighbor_alltoall(sendbuf, sendcount, sendtype, recvbuf,
            recvcount, recvtype, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm, sendcount, sendtype, Neighbor_alltoall,
            t_elapsed, 0);
    } else {
        ret = PMPI_Neighbor_alltoall(sendbuf, sendcount, sendtype, recvbuf,
            recvcount, recvtype, comm);
    }
    return ret;
}

extern "C" {
void mpi_neighbor_alltoall_(const void  *sendbuf, int  *sendcount,
    MPI_Fint  *sendtype, void  *recvbuf, int  *recvcount, MPI_Fint  *recvtype,
                                    MPI_Fint  *comm, MPI_Fint *ierr) {
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    // IN PLACE not meaningful

    ret = MPI_Neighbor_alltoall(sendbuf, *sendcount, c_sendtype, recvbuf,
        *recvcount, c_recvtype, c_comm);
    *ierr = ret;
}
}

int MPI_Neighbor_alltoallv(const void *sendbuf, const int sendcounts[],
    const int sdispls[], MPI_Datatype sendtype, void *recvbuf,
    const int recvcounts[], const int rdispls[], MPI_Datatype recvtype,
    MPI_Comm comm) {
    int ret, sz, i, status, rank, tmp, temp, ndims, j;
    int64_t sum = 0;
    double t_elapsed;
    if ( prof_enabled == 1 ) {
        t_elapsed = MPI_Wtime();
        ret = PMPI_Neighbor_alltoallv(sendbuf, sendcounts, sdispls, sendtype,
            recvbuf, recvcounts, rdispls, recvtype, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;
        // Need to know the topology and the number of neighbors
        PMPI_Comm_rank(comm, &rank);
        PMPI_Topo_test(comm, &status);
        switch (status) {
            case MPI_GRAPH: {
                MPI_Graph_neighbors_count(comm, rank, &sz);
                for ( i=0; i < sz; ++i ) {
                    if ( sendcounts[i] > 0 )
                        sum+=sendcounts[i];
                }
                break;
            }
            case MPI_DIST_GRAPH: {
                MPI_Dist_graph_neighbors_count(comm, &tmp, &sz, &temp);
                for ( i=0; i < sz; ++i ) {
                    if ( sendcounts[i] > 0 )
                        sum+=sendcounts[i];
                }
                break;
            }
            case MPI_CART: {
                // Find the dimensions to determine the max number neighbors
                j = 0;
                PMPI_Cartdim_get(comm, &ndims);
                for ( i = 0; i < ndims; ++i ) {
                    PMPI_Cart_shift(comm, i, 1, &temp, &tmp);
                    if (temp != MPI_PROC_NULL) {
                        sum += sendcounts[j];
                    }
                    j++;
                    if (tmp != MPI_PROC_NULL) {
                        sum += sendcounts[j];
                    }
                    j++;
                }
                break;
            }
            case MPI_KEYVAL_INVALID: {
                mcpt_abort("MPI_TOPO_TEST returned MPI_KEYVAL_INVALID\n");
                break;
            }
            default: {
                mcpt_abort("Unknown Communicator type\n");
                break;
            }
        }
        profile_this(comm, sum, sendtype, Neighbor_alltoallv, t_elapsed, 1);
    } else {
        ret = PMPI_Neighbor_alltoallv(sendbuf, sendcounts, sdispls, sendtype,
            recvbuf, recvcounts, rdispls, recvtype, comm);
    }
    return ret;
}

extern "C" {
void mpi_neighbor_alltoallv_(const void  *sendbuf, const int *sendcounts,
    const int *sdispls, MPI_Fint *sendtype, void *recvbuf,
    const int *recvcounts, const int *rdispls, MPI_Fint *recvtype,
    MPI_Fint *comm, MPI_Fint *ierr) {
    int ret;
    int64_t sum = 0;
    double t_elapsed;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);
    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);


    if ( prof_enabled == 1 ) {
        t_elapsed = MPI_Wtime();
        ret = PMPI_Neighbor_alltoallv(sendbuf, sendcounts, sdispls,
            c_sendtype, recvbuf, recvcounts, rdispls, c_recvtype, c_comm);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(c_comm, sum, c_sendtype, Neighbor_alltoallv,
            t_elapsed, 1);
    } else {
        ret = MPI_Neighbor_alltoallv(sendbuf, sendcounts, sdispls, c_sendtype,
            recvbuf, recvcounts, rdispls, c_recvtype, c_comm);
    }
    *ierr = ret;
}
}

int MPI_Neighbor_alltoallw(const void *sendbuf, const int *sendcounts,
                       const MPI_Aint *sdispls, const MPI_Datatype *sendtypes,
                       void *recvbuf, const int *recvcounts,
                       const MPI_Aint *rdispls, const MPI_Datatype *recvtypes,
                       MPI_Comm comm) {
    int ret, sz, i, j, rank, status, tmp, temp, ndims;
    double t_elapsed;
    int64_t sum = 0;
    if ( prof_enabled == 1 ) {
        t_elapsed = MPI_Wtime();
        ret = PMPI_Neighbor_alltoallw(sendbuf, sendcounts, sdispls, sendtypes,
            recvbuf, recvcounts, rdispls, recvtypes, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;
        PMPI_Comm_rank(comm, &rank);
        PMPI_Topo_test(comm, &status);
        switch (status) {
            case MPI_GRAPH: {
                MPI_Graph_neighbors_count(comm, rank, &sz);
                for ( i = 0; i < sz; i++ ) {
                    if ( sendcounts[i] > 0 ) {
                        PMPI_Type_size(sendtypes[i], &sz);
                        sum+=sendcounts[i]*sz;
                    }
                }
                break;
            }
            case MPI_DIST_GRAPH: {
                MPI_Dist_graph_neighbors_count(comm, &tmp, &sz, &temp);
                for ( i = 0; i < sz; i++ ) {
                    if ( sendcounts[i] > 0 ) {
                        PMPI_Type_size(sendtypes[i], &sz);
                        sum+=sendcounts[i]*sz;
                    }
                }
                break;
            }
            case MPI_CART: {
                // Find the dimensions to determine the max number neighbors
                j = 0;
                PMPI_Cartdim_get(comm, &ndims);
                for ( i =0; i < ndims; ++i ) {
                    PMPI_Cart_shift(comm, i, 1, &temp, &tmp);
                    if (temp != MPI_PROC_NULL) {
                        PMPI_Type_size(sendtypes[j], &sz);
                        sum += sendcounts[j]*sz;
                    }
                    j++;
                    if (tmp != MPI_PROC_NULL) {
                        PMPI_Type_size(sendtypes[j], &sz);
                        sum += sendcounts[j]*sz;
                    }
                    j++;
                }
                break;
            }
            case MPI_KEYVAL_INVALID: {
                mcpt_abort("MPI_TOPO_TEST returned MPI_KEYVAL_INVALID\n");
                break;
            }
            default: {
                mcpt_abort("Unknown Communicator type\n");
                break;
            }
        }
        profile_this(comm, sum, MPI_DATATYPE_NULL, Neighbor_alltoallw,
            t_elapsed, 1);
    } else {
        ret = PMPI_Neighbor_alltoallw(sendbuf, sendcounts, sdispls, sendtypes,
            recvbuf, recvcounts, rdispls, recvtypes, comm);
    }
    return ret;
}

extern "C" {
void mpi_neighbor_alltoallw_(const void *sendbuf, const int *sendcounts,
                             const MPI_Aint *sdispls, const MPI_Fint *sendtypes,
                             void *recvbuf, const int *recvcounts,
                             const MPI_Aint *rdispls, const MPI_Fint *recvtypes,
                             MPI_Fint *comm, MPI_Fint *ierr) {
    int ret, sz, i, j, rank, status, tmp, temp, ndims;
    double t_elapsed;
    int64_t sum = 0;
    MPI_Comm c_comm;
    MPI_Datatype *c_sendtypes, *c_recvtypes;

    c_comm = MPI_Comm_f2c(*comm);
    PMPI_Comm_rank(c_comm, &rank);
    PMPI_Topo_test(c_comm, &status);

    switch (status) {
        // Need to know the topology and the number of neighbors
        case MPI_GRAPH: {
            MPI_Graph_neighbors_count(c_comm, rank, &sz);
            c_sendtypes = static_cast<MPI_Datatype*>(
                malloc(sizeof(MPI_Datatype) * sz));
            c_recvtypes = static_cast<MPI_Datatype*>(
                malloc(sizeof(MPI_Datatype) * sz));
            for ( i=0; i < sz; ++i ) {
                if ( sendcounts[i] > 0 ) {
                    sum+=sendcounts[i];
                    c_sendtypes[i] = MPI_Type_f2c(sendtypes[i]);
                    c_recvtypes[i] = MPI_Type_f2c(recvtypes[i]);
                }
            break;
            }
        }
        case MPI_DIST_GRAPH: {
            MPI_Dist_graph_neighbors_count(c_comm, &tmp, &sz, &temp);
            c_sendtypes = static_cast<MPI_Datatype*>(
                malloc(sizeof(MPI_Datatype) * sz));
            c_recvtypes = static_cast<MPI_Datatype*>(
                malloc(sizeof(MPI_Datatype) * sz));
            for ( i = 0; i < sz; ++i ) {
                if ( sendcounts[i] > 0 ) {
                    sum+=sendcounts[i];
                    c_sendtypes[i] = MPI_Type_f2c(sendtypes[i]);
                    c_recvtypes[i] = MPI_Type_f2c(recvtypes[i]);
                }
            }
            break;
        }
        case MPI_CART: {
            // Find the dimensions to determine the max number neighbors
            j = 0;
            PMPI_Cartdim_get(c_comm, &ndims);
            c_sendtypes = static_cast<MPI_Datatype*>(
                malloc(sizeof(MPI_Datatype) * ndims * 2));
            c_recvtypes = static_cast<MPI_Datatype*>(
                malloc(sizeof(MPI_Datatype) * ndims * 2));
            for ( i = 0; i < ndims; ++i ) {
                PMPI_Cart_shift(c_comm, i, 1, &temp, &tmp);
                if (temp != MPI_PROC_NULL) {
                    sum += sendcounts[j];
                    c_sendtypes[j] = MPI_Type_f2c(sendtypes[j]);
                    c_recvtypes[j] = MPI_Type_f2c(recvtypes[j]);
                }
                j++;
                if (tmp != MPI_PROC_NULL) {
                    sum += sendcounts[j];
                    c_sendtypes[j] = MPI_Type_f2c(sendtypes[j]);
                    c_recvtypes[j] = MPI_Type_f2c(recvtypes[j]);
                }
                j++;
            }
            break;
        }
        case MPI_KEYVAL_INVALID: {
            c_sendtypes = NULL;
            c_recvtypes = NULL;
            mcpt_abort("MPI_TOPO_TEST returned MPI_KEYVAL_INVALID\n");
            break;
        }
        default: {
            c_sendtypes = NULL;
            c_recvtypes = NULL;
            mcpt_abort("Unknown Communicator type\n");
            break;
        }
    }

    if (prof_enabled == 1) {
        t_elapsed = MPI_Wtime();
        ret = MPI_Neighbor_alltoallw(sendbuf, sendcounts, sdispls,
            c_sendtypes, recvbuf, recvcounts, rdispls, c_recvtypes, c_comm);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(c_comm, sum, MPI_DATATYPE_NULL, Neighbor_alltoallw,
            t_elapsed, 1);
    } else {
        ret = MPI_Neighbor_alltoallw(sendbuf, sendcounts, sdispls, c_sendtypes,
            recvbuf, recvcounts, rdispls, c_recvtypes, c_comm);
    }
    free(c_sendtypes);
    free(c_recvtypes);
    *ierr = ret;
}
}

int MPI_Ineighbor_allgather(const void *sendbuf, int sendcount,
                            MPI_Datatype sendtype, void *recvbuf,
                            int recvcount, MPI_Datatype recvtype,
                            MPI_Comm comm, MPI_Request *request) {
    int ret;
    double t_elapsed;
    if ( prof_enabled == 1 ) {
        t_elapsed = MPI_Wtime();
        ret = PMPI_Ineighbor_allgather(sendbuf, sendcount, sendtype, recvbuf,
                                       recvcount, recvtype, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm, sendcount, sendtype, Ineighbor_allgather,
            t_elapsed, 0);
    } else {
        ret = PMPI_Ineighbor_allgather(sendbuf, sendcount, sendtype, recvbuf,
                                       recvcount, recvtype, comm, request);
    }
    return ret;
}

extern "C" {
void mpi_ineighbor_allgather_(const void *sendbuf, int sendcount,
                              MPI_Fint *sendtype, void *recvbuf,
                              int  recvcount, MPI_Fint *recvtype,
                              MPI_Fint *comm, MPI_Fint *request,
                              MPI_Fint *ierr) {
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    ret = MPI_Ineighbor_allgather(sendbuf, sendcount, c_sendtype, recvbuf,
                                  recvcount, c_recvtype, c_comm, &c_request);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
}
}


int MPI_Ineighbor_allgatherv(const void *sendbuf, int sendcount,
    MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
    const int displs[], MPI_Datatype recvtype, MPI_Comm comm,
    MPI_Request *request) {
    int ret;
    double t_elapsed;
    if ( prof_enabled == 1 ) {
        t_elapsed = MPI_Wtime();
        ret = PMPI_Ineighbor_allgatherv(sendbuf, sendcount, sendtype, recvbuf,
            recvcounts, displs, recvtype, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm, sendcount, sendtype, Ineighbor_allgatherv,
            t_elapsed, 1);
    } else {
        ret = PMPI_Ineighbor_allgatherv(sendbuf, sendcount, sendtype, recvbuf,
            recvcounts, displs, recvtype, comm, request);
    }
    return ret;
}

extern "C" {
void mpi_ineighbor_allgatherv_(const void *sendbuf, int sendcount,
                              MPI_Fint *sendtype, void *recvbuf,
                              const int  *recvcounts, const int  *displs,
                              MPI_Fint *recvtype, MPI_Fint *comm,
                              MPI_Fint *request, MPI_Fint *ierr) {
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    ret = MPI_Ineighbor_allgatherv(sendbuf, sendcount, c_sendtype, recvbuf,
        recvcounts, displs, c_recvtype, c_comm, &c_request);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
}
}


int MPI_Ineighbor_alltoall(const void *sendbuf, int sendcount,
    MPI_Datatype sendtype, void *recvbuf, int recvcount,
    MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request) {
    int ret;
    double t_elapsed;
    if ( prof_enabled == 1 ) {
        t_elapsed = MPI_Wtime();
        ret = PMPI_Ineighbor_alltoall(sendbuf, sendcount, sendtype, recvbuf,
                                      recvcount, recvtype, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm, sendcount, sendtype, Ineighbor_alltoall,
            t_elapsed, 0);
    } else {
        ret = PMPI_Ineighbor_alltoall(sendbuf, sendcount, sendtype, recvbuf,
                                      recvcount, recvtype, comm, request);
    }
    return ret;
}

extern "C" {
void mpi_ineighbor_alltoall_(const void *sendbuf, int sendcount,
    MPI_Fint *sendtype, void *recvbuf, int  recvcount, MPI_Fint *recvtype,
    MPI_Fint *comm, MPI_Fint *request, MPI_Fint *ierr) {
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    ret = MPI_Ineighbor_alltoall(sendbuf, sendcount, c_sendtype, recvbuf,
                                  recvcount, c_recvtype, c_comm, &c_request);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
}
}


int MPI_Ineighbor_alltoallv(const void *sendbuf, const int sendcounts[],
                            const int sdispls[], MPI_Datatype sendtype,
                            void *recvbuf, const int recvcounts[],
                            const int rdispls[], MPI_Datatype recvtype,
                            MPI_Comm comm, MPI_Request *request) {
    int ret, rank, sz, i, status, tmp, temp, ndims, j;
    double t_elapsed;
    int64_t sum = 0;
    if ( prof_enabled == 1 ) {
        t_elapsed = MPI_Wtime();
        ret = PMPI_Ineighbor_alltoallv(sendbuf, sendcounts, sdispls, sendtype,
                                       recvbuf, recvcounts, rdispls, recvtype,
                                       comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;
        PMPI_Comm_rank(comm, &rank);
        PMPI_Topo_test(comm, &status);
        switch (status) {
            case MPI_GRAPH: {
                MPI_Graph_neighbors_count(comm, rank, &sz);
                for ( i = 0; i < sz; i++ ) {
                    if ( sendcounts[i] > 0 )
                        sum+=sendcounts[i];
                }
                break;
            }
            case MPI_DIST_GRAPH: {
                MPI_Dist_graph_neighbors_count(comm, &tmp, &sz, &temp);
                for ( i = 0; i < sz; i++ ) {
                    if ( sendcounts[i] > 0 )
                        sum+=sendcounts[i];
                }
                break;
            }
            case MPI_CART: {
                // Find the dimensions to determine the max number neighbors
                j = 0;
                PMPI_Cartdim_get(comm, &ndims);
                for ( i = 0; i < ndims; ++i ) {
                    PMPI_Cart_shift(comm, i, 1, &temp, &tmp);
                    if (temp != MPI_PROC_NULL) {
                        sum += sendcounts[j];
                    }
                    j++;
                    if (tmp != MPI_PROC_NULL) {
                        sum += sendcounts[j];
                    }
                    j++;
                }
                break;
            }
            case MPI_KEYVAL_INVALID: {
                mcpt_abort("MPI_TOPO_TEST returned MPI_KEYVAL_INVALID\n");
                break;
            }
            default: {
                mcpt_abort("Unknown Communicator type\n");
                break;
            }
        }
        profile_this(comm, sum, sendtype, Ineighbor_alltoallv, t_elapsed,
            1);
    } else {
        ret = PMPI_Ineighbor_alltoallv(sendbuf, sendcounts, sdispls, sendtype,
                                       recvbuf, recvcounts, rdispls, recvtype,
                                       comm, request);
    }
    return ret;
}

extern "C" {
void mpi_ineighbor_alltoallv_(const void *sendbuf, const int *sendcounts,
    const int *sdispls, MPI_Fint *sendtype, void *recvbuf,
    const int *recvcounts, const int *rdispls, MPI_Fint *recvtype,
    MPI_Fint *comm, MPI_Fint *request, MPI_Fint *ierr) {
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Ineighbor_alltoallv(sendbuf, sendcounts, sdispls, c_sendtype,
        recvbuf, recvcounts, rdispls, c_recvtype, c_comm, &c_request);

    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
}
}


int MPI_Ineighbor_alltoallw(const void *sendbuf, const int *sendcounts,
                        const MPI_Aint *sdispls, const MPI_Datatype *sendtypes,
                        void *recvbuf, const int *recvcounts,
                        const MPI_Aint *rdispls, const MPI_Datatype *recvtypes,
                        MPI_Comm comm, MPI_Request *request) {
    int ret, sz, i;
    double t_elapsed;
    int64_t sum = 0;
    int status, rank, tmp, temp, ndims, j;
    if ( prof_enabled == 1 ) {
        t_elapsed = MPI_Wtime();
        ret = PMPI_Ineighbor_alltoallw(sendbuf, sendcounts, sdispls,
            sendtypes, recvbuf, recvcounts, rdispls, recvtypes, comm,
            request);
        t_elapsed = MPI_Wtime() - t_elapsed;
        // This not correct for the neighbor alltoallw
        PMPI_Comm_rank(comm, &rank);
        PMPI_Topo_test(comm, &status);
        switch (status) {
            case MPI_GRAPH: {
                MPI_Graph_neighbors_count(comm, rank, &sz);
                for ( i = 0; i < sz; i++ ) {
                    if ( sendcounts[i] > 0 ) {
                        PMPI_Type_size(sendtypes[i], &sz);
                        sum+=sendcounts[i]*sz;
                    }
                }
                break;
            }
            case MPI_DIST_GRAPH: {
                MPI_Dist_graph_neighbors_count(comm, &tmp, &sz, &temp);
                for ( i = 0; i < sz; i++ ) {
                    if ( sendcounts[i] > 0 ) {
                        PMPI_Type_size(sendtypes[i], &sz);
                        sum+=sendcounts[i]*sz;
                    }
                }
                break;
            }
            case MPI_CART: {
                // Find the dimensions to determine the max number neighbors
                j = 0;
                PMPI_Cartdim_get(comm, &ndims);
                for ( i = 0; i < ndims; ++i ) {
                    PMPI_Cart_shift(comm, i, 1, &temp, &tmp);
                    if (temp != MPI_PROC_NULL) {
                        PMPI_Type_size(sendtypes[j], &sz);
                        sum += sendcounts[j]*sz;
                    }
                    j++;
                    if (tmp != MPI_PROC_NULL) {
                        PMPI_Type_size(sendtypes[j], &sz);
                        sum += sendcounts[j]*sz;
                    }
                    j++;
                }
                break;
            }
            case MPI_KEYVAL_INVALID: {
                mcpt_abort("MPI_TOPO_TEST returned MPI_KEYVAL_INVALID\n");
                break;
            }
            default: {
                mcpt_abort("Unknown Communicator type\n");
                break;
            }
        }
        profile_this(comm, sum, MPI_DATATYPE_NULL, Ineighbor_alltoallw,
            t_elapsed, 1);
    } else {
        ret = PMPI_Ineighbor_alltoallw(sendbuf, sendcounts, sdispls, sendtypes,
                                       recvbuf, recvcounts, rdispls, recvtypes,
                                       comm, request);
    }
    return ret;
}

extern "C" {
void mpi_ineighbor_alltoallw_(const void *sendbuf, const int *sendcounts,
    const MPI_Aint *sdispls, const MPI_Fint *sendtypes, void *recvbuf,
    const int *recvcounts, const MPI_Aint *rdispls, const MPI_Fint *recvtypes,
    MPI_Fint *comm, MPI_Fint *request, MPI_Fint *ierr) {
    int ret, sz, i;
    double t_elapsed;
    int64_t sum = 0;
    int status, rank, tmp, temp, ndims, j;
    MPI_Comm c_comm;
    MPI_Datatype *c_sendtypes, *c_recvtypes;
    MPI_Request c_request;

    c_comm = MPI_Comm_f2c(*comm);
    PMPI_Comm_rank(c_comm, &rank);
    PMPI_Topo_test(c_comm, &status);

    switch (status) {
        // Need to know the topology and the number of neighbors
        case MPI_GRAPH: {
            MPI_Graph_neighbors_count(c_comm, rank, &sz);
            c_sendtypes = static_cast<MPI_Datatype*>(
                malloc(sizeof(MPI_Datatype) * sz));
            c_recvtypes = static_cast<MPI_Datatype*>(
                malloc(sizeof(MPI_Datatype) * sz));
            for ( i = 0; i < sz; ++i ) {
                if ( sendcounts[i] > 0 ) {
                    sum+=sendcounts[i];
                    c_sendtypes[i] = MPI_Type_f2c(sendtypes[i]);
                    c_recvtypes[i] = MPI_Type_f2c(recvtypes[i]);
                }
            break;
            }
        }
        case MPI_DIST_GRAPH: {
            MPI_Dist_graph_neighbors_count(c_comm, &tmp, &sz, &temp);
            c_sendtypes = static_cast<MPI_Datatype*>(
                malloc(sizeof(MPI_Datatype) * sz));
            c_recvtypes = static_cast<MPI_Datatype*>(
                malloc(sizeof(MPI_Datatype) * sz));
            for ( i = 0; i < sz; i++ ) {
                if ( sendcounts[i] > 0 ) {
                    sum+=sendcounts[i];
                    c_sendtypes[i] = MPI_Type_f2c(sendtypes[i]);
                    c_recvtypes[i] = MPI_Type_f2c(recvtypes[i]);
                }
            }
            break;
        }
        case MPI_CART: {
            // Find the dimensions to determine the max number neighbors
            j = 0;
            PMPI_Cartdim_get(c_comm, &ndims);
            c_sendtypes = static_cast<MPI_Datatype*>(
                malloc(sizeof(MPI_Datatype) * ndims * 2));
            c_recvtypes = static_cast<MPI_Datatype*>(
                malloc(sizeof(MPI_Datatype) * ndims * 2));
            for ( i = 0; i < ndims; ++i ) {
                PMPI_Cart_shift(c_comm, i, 1, &temp, &tmp);
                if (temp != MPI_PROC_NULL) {
                    sum += sendcounts[j];
                    c_sendtypes[j] = MPI_Type_f2c(sendtypes[j]);
                    c_recvtypes[j] = MPI_Type_f2c(recvtypes[j]);
                }
                j++;
                if (tmp != MPI_PROC_NULL) {
                    sum += sendcounts[j];
                    c_sendtypes[j] = MPI_Type_f2c(sendtypes[j]);
                    c_recvtypes[j] = MPI_Type_f2c(recvtypes[j]);
                }
                j++;
            }
            break;
        }
        case MPI_KEYVAL_INVALID: {
            c_sendtypes = NULL;
            c_recvtypes = NULL;
            mcpt_abort("MPI_TOPO_TEST returned MPI_KEYVAL_INVALID\n");
            break;
        }
        default: {
            c_sendtypes = NULL;
            c_recvtypes = NULL;
            mcpt_abort("Unknown Communicator type\n");
            break;
        }
    }
    if (prof_enabled == 1) {
        t_elapsed = MPI_Wtime();
        ret = PMPI_Ineighbor_alltoallw(sendbuf, sendcounts, sdispls,
            c_sendtypes, recvbuf, recvcounts, rdispls, c_recvtypes,
            c_comm, &c_request);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(c_comm, sum, MPI_DATATYPE_NULL, Ineighbor_alltoallw,
            t_elapsed, 1);
    } else {
        ret = PMPI_Ineighbor_alltoallw(sendbuf, sendcounts, sdispls,
            c_sendtypes, recvbuf, recvcounts, rdispls, c_recvtypes,
            c_comm, &c_request);
    }

    free(c_sendtypes);
    free(c_recvtypes);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
}
}
