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
