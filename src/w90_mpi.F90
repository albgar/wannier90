module w90_mpi

  character(len=1) :: dummy_char
  
#ifdef MPI
  include 'mpif.h'
  
  integer, protected   :: mpi_comm_w90 = mpi_comm_world
  !! Communicator for Wannier90 (settable)

CONTAINS

  subroutine set_w90_comm(comm)
    integer, intent(in) :: comm

    mpi_comm_w90 = comm
  end subroutine set_w90_comm

  subroutine get_w90_comm(comm)
    integer, intent(out) :: comm

    comm = mpi_comm_w90
  end subroutine get_w90_comm
#endif
end module w90_mpi
  

