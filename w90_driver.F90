program w90_driver
  use wannier_m

     call wannier_newlib("gaas",dryrun_in=.true.)
     print *, "Done gaas dryrun"
     call wannier_newlib("gaas")
     print *, "Done gaas"
     call wannier_newlib("lead",post_proc_flag_in=.true.)
     print *, "Done lead preprocessing"
     call wannier_newlib("gaas")
     print *, "Done gaas again"
     
end program w90_driver
