program test_NScool_tov
    use constants_def, only: dp
    use constants_lib
    use nucchem_def
    use nucchem_lib
    use superfluid_lib
    use dStar_eos_lib
    use dStar_crust_def
    use dStar_crust_lib
    use NScool_crust_tov
    
    real(dp) :: lgPstart, lgPend, Mcore, Rcore, Tref
    real(dp), dimension(:), pointer :: y
    integer :: eos_handle, ierr
    
    call constants_init('',ierr)
    call check_okay('constants_init',ierr)
    
    call nucchem_init('../../data',ierr)
    call check_okay('nucchem_init',ierr)
	
    call sf_startup('../../data',ierr)
    call check_okay('sf_startup',ierr)
	
    call sf_load_gaps('ns','gc','t72',ierr)
    call check_okay('sf_load_gaps',ierr)
	
    call dStar_eos_startup('../../data')
    call check_okay('dStar_eos_startup',ierr)
	
    eos_handle = alloc_dStar_eos_handle(ierr)
    call check_okay('alloc_dStar_eos_handle',ierr)
    
    ! switch off the warnings about quantum effects
    call dStar_eos_set_controls(eos_handle,suppress_warnings=.TRUE.)
    
    call dStar_crust_startup('../../data',ierr)
    call check_okay('dStar_atm_startup',ierr)
    
    Tref = 1.0d8
    call dStar_crust_load_table('net_rp',eos_handle, Tref,ierr)
    call check_okay('dStar_crust_load_table',ierr)
    print *,'crust table loaded'
    
    Mcore = 1.6     ! Msun
    Rcore = 11.0    ! km
    lgPstart = 32.5
    lgPend = 27.0
    allocate(y(num_tov_variables))
    print *, 'integrating'
    call tov_integrate(lgPstart, lgPend, Mcore, Rcore, 0.05_dp, y, ierr)
    print *,'done'
    print *,y
    
    call tov_write_crust
    deallocate(y)
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

end program test_NScool_tov
