program test_skyrme
	use constants_def, only: dp
	use constants_lib
	use dStar_core_def, only: rho_saturation
	use brown_skyrme

	integer :: i, ierr
	real(dp) :: rrs(3) = [ 1.0, 3.0, 5.0 ]
	real(dp) :: rho, x, eps, P, muhat, cs2
	call constants_init('',ierr)
	call check_okay('constants_init',ierr)
	call initialize_brown_skyrme(ierr)
	call check_okay('initialize_brown_skyrme',ierr)
! 	do i = skyrme_matter,skyrme_sym
! 		print *,'model: ',skyrme_eos(i)% model
! 		print *,'gamma = ',skyrme_eos(i)% gamma
! 		print *,'    a = ',skyrme_eos(i)% a
! 		print *,'    b = ',skyrme_eos(i)% b
! 		print *,'    c = ',skyrme_eos(i)% c
! 		print *,'    d = ',skyrme_eos(i)% d
! 		print *,'    e = ',skyrme_eos(i)% e
! 		print *,''
! 	end do
	
	x = 0.0
	call do_one
	x = 0.5
	call do_one
contains
	subroutine do_one()
		do i = 1, 3
			rho = rrs(i)*rho_saturation
			call get_skyrme_eos(rho,x,eps,P,muhat,cs2,ierr)
			call check_okay('get_skyrme_eos',ierr)
			print *, rrs(i),eps,P,muhat,cs2
		end do
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
