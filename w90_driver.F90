program w90_driver
  use wannier_m

     call wannier_newlib("gaas")
     print *, "Done gaas"
     call wannier_newlib("lead")
     print *, "Done lead"
     
end program w90_driver
