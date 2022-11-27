!-*- mode: F90 -*-!
!------------------------------------------------------------!
!                                                            !
!                       WANNIER90                            !
!                                                            !
!          The Maximally-Localised Generalised               !
!                 Wannier Functions Code                     !
!                                                            !
! Please cite                                                !
!                                                            !
!  [ref] "Wannier90 as a community code:                     !
!        new features and applications",                     !
!        G. Pizzi et al.,  J. Phys. Cond. Matt. 32,          !
!        165902 (2020).                                      !
!        http://doi.org/10.1088/1361-648X/ab51ff             !
!                                                            !
! in any publications arising from the use of this code.     !
!                                                            !
! Wannier90 is based on Wannier77, written by N. Marzari,    !
! I. Souza and D. Vanderbilt. For the method please cite     !
!                                                            !
! [ref] N. Marzari and D. Vanderbilt,                        !
!       Phys. Rev. B 56 12847 (1997)                         !
!       http://dx.doi.org/10.1103/PhysRevB.56.12847          !
!                                                            !
! [ref] I. Souza, N. Marzari and D. Vanderbilt,              !
!       Phys. Rev. B 65 035109 (2001)                        !
!       http://dx.doi.org/10.1103/PhysRevB.65.035109         !
!                                                            !
! [ref] N. Marzari, A. A. Mostofi, J. R. Yates, I. Souza,    !
!       D. Vanderbilt, "Maximally localized Wannier          !
!       functions: theory and applications",                 !
!       Rev. Mod. Phys. 84, 1419 (2012)                      !
!       http://dx.doi.org/10.1103/RevModPhys.84.1419         !
!                                                            !
! For a full list of authors and contributors, please        !
! see the README file in the root directory of the           !
! distribution.                                              !
!                                                            !
! This file is distributed as part of the Wannier90 code and !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the Wannier90      !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the Wannier90 code is www.wannier.org       !
!                                                            !
! The Wannier90 code is hosted on GitHub:                    !
!                                                            !
! https://github.com/wannier-developers/wannier90            !
!------------------------------------------------------------!

module wannier_m

  public :: wannier_newlib
  private
  
CONTAINS
  
  subroutine wannier_newlib(seedname_in, nnkp_mode, dryrun_mode,   &
                            eigfile_ext,                           &
                            nntot_out, nnlist_out, nncell_out,     &
                            u_matrix_out, u_matrix_opt_out)
  
  !! The main Wannier90 program, wrapped as a subroutine

  use w90_constants
  use w90_parameters
  use w90_io
  use w90_hamiltonian
  use w90_kmesh
  use w90_disentangle
  use w90_overlap
  use w90_wannierise
  use w90_plot
  use w90_transport
  use w90_comms, only: on_root, num_nodes, comms_setup, comms_end, comms_bcast, my_node_id
  use w90_sitesym !YN:

  implicit none

  character(len=*), intent(in)   :: seedname_in
  logical, intent(in), optional  :: nnkp_mode
  logical, intent(in), optional  :: dryrun_mode
  character(len=*), intent(in), optional :: eigfile_ext

  ! Information on k-point neighbors that clients can request
  integer, intent(out), optional :: nntot_out
  integer, intent(out), allocatable, optional :: nnlist_out(:,:)
  integer, intent(out), allocatable, optional :: nncell_out(:,:,:)

  ! Information on unitary matrices that clients can request
  complex(dp), intent(out), allocatable, optional :: u_matrix_out(:,:,:)
  complex(dp), intent(out), allocatable, optional :: u_matrix_opt_out(:,:,:)

  real(kind=dp) time0, time1, time2
  character(len=9) :: stat, pos, cdate, ctime
  logical :: wout_found, dryrun
  integer :: len_seedname
  character(len=50) :: prog

  call comms_setup

  library = .false.

  time0 = io_time()

  if (on_root) then
    prog = 'wannier90-newlib'
    !    call io_commandline(prog, dryrun)

    seedname = seedname_in

    if (present(nnkp_mode)) then
       post_proc_flag = nnkp_mode
    else
       post_proc_flag = .false.
    endif
    if (present(dryrun_mode)) then
       dryrun = dryrun_mode
    else
       dryrun = .false.
    endif

    if (present(eigfile_ext)) then
       eig_ext = eigfile_ext
    else
       eig_ext = ".eig"
    endif

    len_seedname = len(seedname)
  end if
 
  call comms_bcast(len_seedname, 1)
  call comms_bcast(seedname, len_seedname)
  call comms_bcast(dryrun, 1)

  if (on_root) then
    stdout = io_file_unit()
    open (unit=stdout, file=trim(seedname)//'.werr')
    call io_date(cdate, ctime)
    write (stdout, *) 'Wannier90: Execution started on ', cdate, ' at ', ctime

    call param_read
    close (stdout, status='delete')

    if (restart .eq. ' ') then
      stat = 'replace'
      pos = 'rewind'
    else
      inquire (file=trim(seedname)//'.wout', exist=wout_found)
      if (wout_found) then
        stat = 'old'
      else
        stat = 'replace'
      endif
      pos = 'append'
    endif

    stdout = io_file_unit()
    open (unit=stdout, file=trim(seedname)//'.wout', status=trim(stat), position=trim(pos))
    call param_write_header()
    if (num_nodes == 1) then
#ifdef MPI
      write (stdout, '(/,1x,a)') 'Running in serial (with parallel executable)'
#else
      write (stdout, '(/,1x,a)') 'Running in serial (with serial executable)'
#endif
    else
      write (stdout, '(/,1x,a,i3,a/)') &
        'Running in parallel on ', num_nodes, ' CPUs'
    endif
    call param_write()

    time1 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') 'Time to read parameters  ', time1 - time0, ' (sec)'

    if (.not. explicit_nnkpts) call kmesh_get
    time2 = io_time()
    write (stdout, '(1x,a25,f11.3,a)') &
      'Time to get kmesh        ', time2 - time1, ' (sec)'

    call param_memory_estimate
  end if

  if (dryrun) then
    if (on_root) then
      write (stdout, *) ' '
      write (stdout, *) '                       ==============================='
      write (stdout, *) '                                   DRYRUN             '
      write (stdout, *) '                       No problems found with win file'
      write (stdout, *) '                       ==============================='
    endif
    call kmesh_dealloc()
    call param_dealloc()
    close(unit=stdout)
    return
  endif

  ! We now distribute the parameters to the other nodes
  call param_dist
  if (gamma_only .and. num_nodes > 1) &
    call io_error('Gamma point branch is serial only at the moment')

  if (transport .and. tran_read_ht) goto 3003

  ! Sort out restarts
  if (restart .eq. ' ') then  ! start a fresh calculation
    if (on_root) write (stdout, '(1x,a/)') 'Starting a new Wannier90 calculation ...'
  else                      ! restart a previous calculation
    if (on_root) call param_read_chkpt()
    call param_chkpt_dist
    if (lsitesymmetry) call sitesym_read()   ! update this to read on root and bcast - JRY

    select case (restart)
    case ('default')    ! continue from where last checkpoint was written
      if (on_root) write (stdout, '(/1x,a)', advance='no') 'Resuming a previous Wannier90 calculation '
      if (checkpoint .eq. 'postdis') then
        if (on_root) write (stdout, '(a/)') 'from wannierisation ...'
        goto 1001         ! go to wann_main
      elseif (checkpoint .eq. 'postwann') then
        if (on_root) write (stdout, '(a/)') 'from plotting ...'
        goto 2002         ! go to plot_main
      else
        if (on_root) write (stdout, '(/a/)')
        call io_error('Value of checkpoint not recognised in wann_prog')
      endif
    case ('wannierise') ! continue from wann_main irrespective of value of last checkpoint
      if (on_root) write (stdout, '(1x,a/)') 'Restarting Wannier90 from wannierisation ...'
      goto 1001
    case ('plot')       ! continue from plot_main irrespective of value of last checkpoint
      if (on_root) write (stdout, '(1x,a/)') 'Restarting Wannier90 from plotting routines ...'
      goto 2002
    case ('transport')   ! continue from tran_main irrespective of value of last checkpoint
      if (on_root) write (stdout, '(1x,a/)') 'Restarting Wannier90 from transport routines ...'
      goto 3003
    case default        ! for completeness... (it is already trapped in param_read)
      call io_error('Value of restart not recognised in wann_prog')
    end select
  endif

  if (postproc_setup) then
    if (on_root) call kmesh_write()
    !
    !
    if (present(nntot_out)) then
       nntot_out = nntot
    endif
    if (present(nnlist_out)) then
       if (allocated(nnlist_out)) deallocate(nnlist_out)
       allocate(nnlist_out,source=nnlist)
    endif
    if (present(nncell_out)) then
       if (allocated(nncell_out)) deallocate(nncell_out)
       allocate(nncell_out,source=nncell)
    endif
    !
    !  sanity check
    !
    if (present(u_matrix_out) .or. present(u_matrix_opt_out)) then
       call io_error('Error: Cannot request u_matrices in nnkp mode')
    endif
    !
    call kmesh_dealloc()
    call param_dealloc()
    if (on_root) write (stdout, '(1x,a25,f11.3,a)') 'Time to write kmesh      ', io_time(), ' (sec)'
    if (on_root) write (stdout, '(/a)') ' Exiting... '//trim(seedname)//'.nnkp written.'
    call comms_end
    return
  endif

  if (lsitesymmetry) call sitesym_read()   ! update this to read on root and bcast - JRY
  call overlap_allocate()
  call overlap_read()

  time1 = io_time()
  if (on_root) write (stdout, '(/1x,a25,f11.3,a)') 'Time to read overlaps    ', time1 - time2, ' (sec)'

  have_disentangled = .false.

  if (disentanglement) then
    call dis_main()
    have_disentangled = .true.
    time2 = io_time()
    if (on_root) write (stdout, '(1x,a25,f11.3,a)') 'Time to disentangle bands', time2 - time1, ' (sec)'
  endif

  if (on_root) call param_write_chkpt('postdis')
!~  call param_write_um

1001 time2 = io_time()

  if (.not. gamma_only) then
    call wann_main()
  else
    call wann_main_gamma()
  end if

  time1 = io_time()
  if (on_root) write (stdout, '(1x,a25,f11.3,a)') 'Time for wannierise      ', time1 - time2, ' (sec)'

  if (on_root) call param_write_chkpt('postwann')

2002 continue
  if (on_root) then
    ! I call the routine always; the if statements to decide if/what
    ! to plot are inside the function
    time2 = io_time()
    call plot_main()
    time1 = io_time()
    ! Now time is always printed, even if no plotting is done/required, but
    ! it shouldn't be a problem.
    write (stdout, '(1x,a25,f11.3,a)') 'Time for plotting        ', time1 - time2, ' (sec)'
  endif

3003 continue
  if (on_root) then
    time2 = io_time()
    if (transport) then
      call tran_main()
      time1 = io_time()
      write (stdout, '(1x,a25,f11.3,a)') 'Time for transport       ', time1 - time2, ' (sec)'
      if (tran_read_ht) goto 4004
    end if
  endif

  ! Transfer info on unitary matrices if requested
  if (present(u_matrix_out)) then
     if (allocated(u_matrix_out)) deallocate(u_matrix_out)
     allocate(u_matrix_out,source=u_matrix)
  endif
  if (present(u_matrix_opt_out)) then
     if (allocated(u_matrix_opt_out)) deallocate(u_matrix_opt_out)
     allocate(u_matrix_opt_out,source=u_matrix_opt)
  endif

  call tran_dealloc()
  call hamiltonian_dealloc()
  call overlap_dealloc()
  call kmesh_dealloc()
  call param_dealloc()
  if (lsitesymmetry) call sitesym_dealloc() !YN:

4004 continue

  if (on_root) then
    write (stdout, '(1x,a25,f11.3,a)') 'Total Execution Time     ', io_time(), ' (sec)'

    if (timing_level > 0) call io_print_timings()

    write (stdout, *)
    write (stdout, '(1x,a)') 'All done: wannier90 exiting'

    close (stdout)
  endif

  call comms_end

end subroutine wannier_newlib

end module wannier_m
