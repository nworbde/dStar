program test_atm
    use constants_lib
    use constants_lib
    use nucchem_lib
    use dStar_eos_lib
    use conductivity_lib
	use dStar_atm_def
	use dStar_atm_lib
	
	real(dp) :: lgTb,g,Plight,Pb,lgTeff,dlgTeff,lgflux,dlgflux
	integer :: ierr,i
    
    call constants_init('',ierr)
	call check_okay('constants_init',ierr)
    call nucchem_init('../../data',ierr)
	call check_okay('nucchem_init',ierr)
    call dStar_eos_startup('../../data')
	call dStar_atm_startup('../../data',ierr)
	call check_okay('dStar_atm_startup',ierr)
	
	write (*,'(5(a9,tr2),/,5("=========",tr2))') 'lgTb','lgTeff','dlgTeff','lgflux','dlgflux'
    
	g = 2.43e14_dp
    Pb = 4.3e13_dp*g
    
	Plight = 2.43e19_dp    
    call do_one('bc09')
    call do_one('pcy97')
    
	Plight = 2.43e23_dp
    call do_one('bc09')
    call do_one('pcy97')    
    
	call dStar_atm_shutdown
    call dStar_eos_shutdown
    call nucchem_shutdown
    
contains
    
    subroutine do_one(prefix)
        character(len=*), intent(in) :: prefix
        integer :: i

        write (*,'(/,a)') trim(prefix)
        write (*,'(a,f6.2)') 'Plight = ',log10(Plight)
        
    	call dStar_atm_load_table(trim(prefix),g,Plight,Pb,ierr)
    	call check_okay('dStar_atm_load_table',ierr)
        
    	do i = 1, 20
    		lgTb = (7.0_dp + 1.5_dp*real(i-1.0,dp)/19.0_dp)
    		call dStar_atm_get_results(lgTb,lgTeff,dlgTeff,lgflux,dlgflux,ierr)
    		call check_okay('dStar_atm_get_results',ierr)
    		write (*,'(5(f9.5,tr2))') lgTb,lgTeff,dlgTeff,lgflux,dlgflux
    	end do
    	call dStar_atm_free_table
        
    end subroutine do_one
    
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
