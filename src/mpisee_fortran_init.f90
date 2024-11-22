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
