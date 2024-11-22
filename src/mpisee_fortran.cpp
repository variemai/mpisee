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
