! fortran_constants.f90
!module fortran_constants
!  implicit none

!contains

  subroutine mpisee_fortran_in_place_init() bind(C)
  use mpi
  call mpisee_init_fortran_mpi_in_place(MPI_IN_PLACE)
  end subroutine mpisee_fortran_in_place_init

!end module fortran_constants
