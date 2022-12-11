
set(test_dir "${CMAKE_BINARY_DIR}/install_tests")

file(MAKE_DIRECTORY "${test_dir}")

function(setup_test name)
 file(COPY "${name}/"  DESTINATION "${test_dir}")
endfunction()
#
setup_test("example01")
setup_test("example02")
setup_test("example03")
setup_test("example04")

if(WITH_MPI)
   set(num_ranks 2)
   set(mpiexec_args
          ${MPIEXEC_NUMPROC_FLAG} ${num_ranks}
          ${MPIEXEC_PREFLAGS}
	  $<TARGET_FILE:w90_driver> 
          ${MPIEXEC_POSTFLAGS})
	  
   add_test( NAME    w90-driver_mpi
             COMMAND ${MPIEXEC_EXECUTABLE} ${mpiexec_args}
             WORKING_DIRECTORY "${test_dir}" )

 else()
 
   add_test( NAME    w90-driver
             COMMAND 	  $<TARGET_FILE:w90_driver> 
             WORKING_DIRECTORY "${test_dir}" )

endif(WITH_MPI)
