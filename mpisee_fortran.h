#ifndef MPISEE_FORTRAN_H_
#define MPISEE_FORTRAN_H_
#include "mpi.h"

extern void *mpisee_fortran_mpi_in_place;

extern "C" void mpisee_init_fortran_mpi_in_place_(MPI_Fint *in_place);
extern "C" void mpisee_fortran_in_place_init( void );

#endif // MPISEE_FORTRAN_H_
