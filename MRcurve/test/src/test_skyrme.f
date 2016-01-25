program test_skyrme
	use constants_def, only: dp
	use constants_lib
	use dStar_core_def
	use dStar_core_lib
	use load_skyrme_parameters
	use brown_skyrme

	integer, parameter :: Ntab = 20
	integer :: i, ierr
	real(dp), dimension(Ntab) :: rho, x, eps, P, muhat, mue, cs2
	
	call constants_init('',ierr)
	call check_okay('constants_init',ierr)
	call dstar_core_startup('../../data',ierr)
	call check_okay('dstar_core_startup',ierr)
	call load_skyrme_table('s7d',ierr)
	call check_okay('load_skyrme_table',ierr)
	
	rho = [(0.5 + 6.0*real(i-1,dp)/real(Ntab-1,dp), i=1,Ntab)]
	rho = rho*rho_saturation
	call find_beta_eq(rho,x,eps,P,muhat,mue,cs2,ierr)
	call check_okay('find_beta_eq',ierr)

	do i = 1,Ntab
		write(*,*) rho(i), x(i), eps(i), P(i), cs2(i), muhat(i)-mue(i)
	end do
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

end program test_skyrme
