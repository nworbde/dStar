program test_TOV
	use constants_def
	use constants_lib
    use nucchem_def
    use nucchem_lib
    use superfluid_lib
    use dStar_eos_lib
    use dStar_crust_def
    use dStar_crust_lib
	use dStar_core_def
	use dStar_core_lib
	use dStar_core_tov

	character(len=*), parameter :: datadir='../../data'
	integer :: ierr, i
    integer :: eos_handle
    real(dp) :: Tref, lgPstart, lgPend, lnP_rez
	real(dp), dimension(:), pointer :: y
	real(dp), dimension(:), allocatable :: M, R
	
	call constants_init('',ierr)
	call check_okay('constants_init',ierr)

	call dstar_core_startup(datadir,ierr)
	call check_okay('dstar_core_startup',ierr)
    
    call nucchem_init(datadir,ierr)
    call check_okay('nucchem_init',ierr)
	
    call sf_startup(datadir,ierr)
    call check_okay('sf_startup',ierr)
	
    call sf_load_gaps('ns','gc','t72',ierr)
    call check_okay('sf_load_gaps',ierr)
	
    call dStar_eos_startup(datadir)
    call check_okay('dStar_eos_startup',ierr)
	
    eos_handle = alloc_dStar_eos_handle(ierr)
    call check_okay('alloc_dStar_eos_handle',ierr)
    
    ! switch off the warnings about quantum effects
    call dStar_eos_set_controls(eos_handle,suppress_warnings=.TRUE.)
    
    call dStar_crust_startup(datadir,ierr)
    call check_okay('dStar_atm_startup',ierr)
    
    Tref = 1.0d8
    call dStar_crust_load_table('hz90',eos_handle, Tref,ierr)
    call check_okay('dStar_crust_load_table',ierr)
	
	call dStar_core_load_table('s7d',ierr)
	call check_okay('dStar_core_load_table',ierr)
	
	lgPend = -7.0
	lnP_rez = 0.2
	allocate(y(num_core_variables))
	allocate(M(20),R(20))
	do i = 1,20
		lgPstart = 1.1 + 2.2*real(i-1,dp)/19.0
		call core_integrate(lgPstart, lgPend, lnP_rez, y, ierr)
		call check_okay('core_integrate',ierr)
		M(i) = y(core_mass)
		R(i) = y(core_radius)*length_g*1.0e-5
		Pc(i) = core_structure% P(core_structure% nzs)
		call core_write_crust
	end do

	write (*,*)
	do i = 1, 20
		write (*,'(2f6.2)') M(i), R(i)
	end do

	call dStar_core_shutdown
    call dStar_crust_shutdown
    call sf_shutdown
    call nucchem_shutdown
	
contains
	subroutine check_okay(msg,ierr)
		use iso_fortran_env, only : error_unit
		character(len=*), intent(in) :: msg
		integer, intent(inout) :: ierr
		if (ierr /= 0) then
			write (error_unit,*) trim(msg)//': ierr = ',ierr
			if (ierr < 0) stop
		end if
	end subroutine check_okay

end program test_TOV
