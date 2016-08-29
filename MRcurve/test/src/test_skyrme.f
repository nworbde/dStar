program test_skyrme
	use constants_def
	use constants_lib
	use dStar_core_def
	use dStar_core_lib

	integer :: i, ierr
	real(dp) :: lgP,lgRho,lgEps,Xprot,dlgRho,dlgEps,dXprot
	
	call constants_init('',ierr)
	call check_okay('constants_init',ierr)
	call dstar_core_startup('../../data',ierr)
	call check_okay('dstar_core_startup',ierr)

	write (*,'(7(a9,tr2),/,7("=========",tr2))') 'lgP','lgRho','dlgRho', &
		& 'lgEps','dlgEps','Xprot','dXprot'

	call do_one('s7d')
	call do_one('s7a')
	call do_one('s7r')
	call dStar_core_shutdown
contains
	subroutine do_one(model)
		character(len=*), intent(in) :: model
		integer :: i
		
        write (*,'(/,a)') trim(model)
        
    	call dStar_core_load_table(trim(model),ierr)
    	call check_okay('dStar_core_load_table',ierr)
        
    	do i = 1, 20
    		lgP = (-0.5_dp + 3.5_dp*real(i-1.0,dp)/19.0_dp)
    		call dStar_core_get_results(lgP,lgRho,dlgRho,lgEps,dlgEps,Xprot,dXprot,ierr)
    		call check_okay('dStar_core_get_results',ierr)
    		write (*,'(7(f9.5,tr2))') lgP,lgRho,dlgRho,lgEps,dlgEps,Xprot,dXprot
    	end do
    	call dStar_core_free_table
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

end program test_skyrme
