program w90_driver
  use wannier_m
#ifdef MPI
  use mpi
#endif
  
  integer, parameter :: dp = selected_real_kind(10)
  
  integer :: nntot
  integer,  allocatable :: nncell(:,:,:), nnlist(:,:)

#ifdef MPI
  integer MPI_err
  call MPI_init(MPI_err)
#endif
  
  call wannier_newlib("minimal",         &
#ifdef MPI
                       mpi_comm=mpi_comm_world, &
#endif
                       nnkp_mode=.true., &
                       nntot_out=nntot, &
                       nnlist_out=nnlist, &
                       nncell_out=nncell)
                       
     print *, "Done minimal preprocessing"
     print *, "nntot: ", nntot
     print *, "shape nnlist: ", shape(nnlist)
     print *, "shape nncell: ", shape(nncell)
     
     call wannier_newlib("gaas",  &
#ifdef MPI
                       mpi_comm=mpi_comm_world, &
#endif
                       dryrun_mode=.true.)
     print *, "Done gaas dryrun"
     call wannier_newlib("gaas", &
#ifdef MPI
                       mpi_comm=mpi_comm_world &
#endif
                        )
     print *, "Done gaas"

     call wannier_newlib("lead",  &
#ifdef MPI
                       mpi_comm=mpi_comm_world, &
#endif
                       nnkp_mode=.true., &               
                       nntot_out=nntot, &
                       nnlist_out=nnlist, &
                       nncell_out=nncell)
     print *, "Done lead preprocessing"
     print *, "nntot: ", nntot
     print *, "shape nnlist: ", shape(nnlist)
     print *, "shape nncell: ", shape(nncell)

     call wannier_newlib("gaas", &
#ifdef MPI
                       mpi_comm=mpi_comm_world &
#endif
                        )
     print *, "Done gaas again"

     call wannier_newlib("../example04/copper", &
#ifdef MPI
                       mpi_comm=mpi_comm_world &
#endif
                        )
     print *, "Done copper"
     
#ifdef MPI
  call MPI_finalize(MPI_err)
#endif
  
end program w90_driver
