program w90_driver
  use wannier_m
  integer, parameter :: dp = selected_real_kind(10)
  
  integer :: nntot
  integer,  allocatable :: nncell(:,:,:), nnlist(:,:)
  
  call wannier_newlib("minimal",nnkp_mode=.true., &
                       nntot_out=nntot, &
                       nnlist_out=nnlist, &
                       nncell_out=nncell)
                       
     print *, "Done minimal preprocessing"
     print *, "nntot: ", nntot
     print *, "shape nnlist: ", shape(nnlist)
     print *, "shape nncell: ", shape(nncell)
     
     call wannier_newlib("gaas",dryrun_mode=.true.)
     print *, "Done gaas dryrun"
     call wannier_newlib("gaas")
     print *, "Done gaas"

     call wannier_newlib("lead",nnkp_mode=.true., &
                       nntot_out=nntot, &
                       nnlist_out=nnlist, &
                       nncell_out=nncell)
     print *, "Done lead preprocessing"
     print *, "nntot: ", nntot
     print *, "shape nnlist: ", shape(nnlist)
     print *, "shape nncell: ", shape(nncell)

     call wannier_newlib("gaas")
     print *, "Done gaas again"
     
end program w90_driver
