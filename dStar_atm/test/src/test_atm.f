program test_atm
    use constants_lib
	use dStar_atm_def
	use dStar_atm_lib
	
	real(dp) :: lgTb,g,Plight,lgTeff,dlgTeff,lgflux,dlgflux
	integer :: ierr,i
    
    call constants_init('',ierr)
	call dStar_atm_startup('../../data',ierr)
	call check_okay('dStar_atm_startup',ierr)
	
	g = 2.43e14_dp
	Plight = 1.0e6_dp
	call dStar_atm_load_table('pcy97',g,Plight,ierr)
	call check_okay('dStar_atm_load_table',ierr)
	do i = 1, 20
		lgTb = (7.0_dp + 1.5_dp*real(i-1.0,dp)/19.0_dp)
		call dStar_atm_get_results(lgTb,lgTeff,dlgTeff,lgflux,dlgflux, ierr)
		call check_okay('dStar_atm_get_results',ierr)
		write (*,'(5(f9.5,tr2))') lgTb,lgTeff,dlgTeff,lgflux,dlgflux
	end do
	
	call dStar_atm_free_table
	Plight = 1.0e28_dp
	call dStar_atm_load_table('pcy97',g,Plight,ierr)
	call check_okay('dStar_atm_load_table',ierr)
	do i = 1, 20
		lgTb = (7.0_dp + 1.5_dp*real(i-1.0,dp)/19.0_dp)
		call dStar_atm_get_results(lgTb,lgTeff,dlgTeff,lgflux,dlgflux, ierr)
		call check_okay('dStar_atm_get_results',ierr)
		write (*,'(5(f9.5,tr2))') lgTb,lgTeff,dlgTeff,lgflux,dlgflux
	end do
	
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
end program test_atm
