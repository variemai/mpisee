#include "mpisee_fortran.h"
#include <stdio.h>

void *mpisee_fortran_mpi_in_place = NULL;
void *mpisee_fortran_mpi_status_ignore = NULL;
void *mpisee_fortran_mpi_statuses_ignore = NULL;
void *mpisee_fortran_mpi_bottom = NULL ;
void *mpisee_fortran_mpi_unweighted = NULL;

extern "C" void mpisee_init_fortran_mpi_in_place_(MPI_Fint *in_place) {
    mpisee_fortran_mpi_in_place = in_place;
    //printf("mpisee: fortran_mpi_in_place is set to %p\n", fortran_mpi_in_place);
}

extern "C" void mpisee_init_fortran_mpi_status_ignore_(MPI_Fint *status_ignore) {
    mpisee_fortran_mpi_status_ignore = status_ignore;
    //printf("mpisee: fortran_mpi_status_ignore is set to %p\n", fortran_mpi_status_ignore);
}

extern "C" void mpisee_init_fortran_mpi_statuses_ignore_(MPI_Fint *statuses_ignore) {
    mpisee_fortran_mpi_statuses_ignore = statuses_ignore;
    //printf("mpisee: fortran_mpi_statuses_ignore is set to %p\n", fortran_mpi_statuses_ignore);
}

extern "C" void mpisee_init_fortran_mpi_bottom_(MPI_Fint *bottom) {
    mpisee_fortran_mpi_bottom = bottom;
    //printf("mpisee: fortran_mpi_bottom is set to %p\n", fortran_mpi_bottom);
}

extern "C" void mpisee_init_fortran_mpi_unweighted_(MPI_Fint *unweighted) {
    mpisee_fortran_mpi_unweighted = unweighted;
    //printf("mpisee: fortran_mpi_unweighted is set to %p\n", fortran_mpi_unweighted);
}
