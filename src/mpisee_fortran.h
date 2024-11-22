#ifndef MPISEE_FORTRAN_H_
#define MPISEE_FORTRAN_H_
#include "mpi.h"

extern void *mpisee_fortran_mpi_in_place;
extern void *mpisee_fortrtan_mpi_status_ignore;
extern void *mpisee_fortrtan_mpi_statuses_ignore;
extern void *mpisee_fortran_mpi_bottom;
extern void *mpisee_fortran_mpi_unweighted;

extern "C" void mpisee_init_fortran_mpi_in_place_(MPI_Fint *in_place);
extern "C" void mpisee_fortran_in_place_init( void );
extern "C" void mpisee_init_fortran_mpi_status_ignore_(MPI_Fint *status_ignore);
extern "C" void mpisee_fortran_status_ignore_init( void );
extern "C" void mpisee_init_fortran_mpi_statuses_ignore_(MPI_Fint *statuses_ignore);
extern "C" void mpisee_fortran_statuses_ignore_init( void );
extern "C" void mpisee_init_fortran_mpi_bottom_(MPI_Fint *bottom);
extern "C" void mpisee_fortran_bottom_init( void );
extern "C" void mpisee_init_fortran_mpi_unweighted_(MPI_Fint *unweighted);
extern "C" void mpisee_fortran_unweighted_init( void );


#endif // MPISEE_FORTRAN_H_
