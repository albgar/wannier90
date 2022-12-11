program w90_driver
  use wannier90_m, only: wannier90_wrapper
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
  
     call wannier90_wrapper("gaas",  &
#ifdef MPI
                       mpi_comm=mpi_comm_world, &
#endif
                       dryrun_mode=.true.)
     print *, "Done gaas dryrun"
     
#ifdef MPI
     call wannier90_wrapper("gaas", mpi_comm=mpi_comm_world)
#else
     call wannier90_wrapper("gaas")
#endif
                        
     print *, "Done gaas"

#ifdef MPI
     call wannier90_wrapper("silicon", mpi_comm=mpi_comm_world)
#else
     call wannier90_wrapper("silicon")
#endif
                        
     print *, "Done silicon"

     call wannier90_wrapper("lead",  &
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

#ifdef MPI
     call wannier90_wrapper("gaas", mpi_comm=mpi_comm_world)
#else
     call wannier90_wrapper("gaas")
#endif

     print *, "Done gaas again"

#ifdef MPI
     call wannier90_wrapper("copper", mpi_comm=mpi_comm_world)
#else
     call wannier90_wrapper("copper")
#endif
     print *, "Done copper"
     
#ifdef MPI
  call MPI_finalize(MPI_err)
#endif
  
end program w90_driver
