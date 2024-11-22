!******************************************************************************
! This source file is part of mpisee
! Copyright (C) 2024  Ioannis Vardas - vardas@par.tuwien.ac.at
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>
!******************************************************************************

! fortran_constants.f90
!module fortran_constants
!  implicit none

!contains

  subroutine mpisee_fortran_in_place_init() bind(C)
  use mpi
  call mpisee_init_fortran_mpi_in_place(MPI_IN_PLACE)
  end subroutine mpisee_fortran_in_place_init

  subroutine mpisee_fortran_status_ignore_init() bind(C)
  use mpi
  call mpisee_init_fortran_mpi_status_ignore(MPI_STATUS_IGNORE)
  end subroutine mpisee_fortran_status_ignore_init

  subroutine mpisee_fortran_statuses_ignore_init() bind(C)
  use mpi
  call mpisee_init_fortran_mpi_statuses_ignore(MPI_STATUSES_IGNORE)
  end subroutine mpisee_fortran_statuses_ignore_init

  subroutine mpisee_fortran_bottom_init() bind(C)
  use mpi
  call mpisee_init_fortran_mpi_bottom(MPI_BOTTOM)
  end subroutine mpisee_fortran_bottom_init

  subroutine mpisee_fortran_unweighted_init() bind(C)
  use mpi
  call mpisee_init_fortran_mpi_unweighted(MPI_UNWEIGHTED)
  end subroutine mpisee_fortran_unweighted_init

!end module fortran_constants
