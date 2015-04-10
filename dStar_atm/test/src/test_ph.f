program test_ph
    use constants_lib
	use dStar_atm_def
	use dStar_atm_lib
	use bc09
	
	real(dp) :: lgTb,g,Plight,lgTeff,dlgTeff,lgflux,dlgflux
	integer :: ierr,i
    
    call constants_init('',ierr)
	call dStar_atm_startup('../../data',ierr)
	call check_okay('dStar_atm_startup',ierr)
	
	g = 2.43e14_dp
	Plight = 1.0e6_dp
	
	call do_get_bc09_Teff(grav, Plight, Tb, Teff, flux)
	call check_okay('do_get_bc09_Teff',ierr)

	call dStar_atm_shutdown
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
end program test_ph
