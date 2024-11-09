#include "mpisee_fortran.h"
#include <stdio.h>

void *mpisee_fortran_mpi_in_place = NULL;

extern "C" void mpisee_init_fortran_mpi_in_place_(MPI_Fint *in_place) {
    mpisee_fortran_mpi_in_place = in_place;
    //printf("mpisee: fortran_mpi_in_place is set to %p\n", fortran_mpi_in_place);
}
