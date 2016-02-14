module beta_equilibrium
	use constants_def
	use dStar_core_def
contains
	
	subroutine find_beta_eq(eos_handle,rho,T,x,eps,P,muhat,mue,ierr)
		use constants_def
		use dStar_eos_def
		use dStar_eos_lib
		use num_lib
		integer, intent(in) :: eos_handle
		real(dp), dimension(:), intent(in) :: rho
		real(dp), intent(in) :: T
		real(dp), dimension(:), intent(out) :: x,eps,P,muhat,mue
		integer, intent(out) :: ierr
		integer, parameter :: number_integer_parameters = 1
		integer, parameter :: number_real_parameters = 6
        real(dp), dimension(:), pointer :: rpar=>null()
        integer, dimension(:), pointer :: ipar=>null()
        integer :: lipar, lrpar
        integer :: i,Ntab
        real(dp) :: x1, x3, y1, y3, epsx, epsy, xguess, dfdx
        integer :: imax
        
        Ntab = size(rho)
        imax = 20
        epsx = 1.0d-8
        epsy = 1.0d-8
		
		lipar = number_integer_parameters
		allocate(ipar(lipar))
		ipar(1) = eos_handle
		
		lrpar = number_real_parameters
		allocate(rpar(lrpar))
		
		do i = 1, Ntab
			rpar(1) = rho(i)
			rpar(2) = T
			x3 = 0.5_dp
			y3 = muhat_mue(x3,dfdx,lrpar,rpar,lipar,ipar,ierr)
			x1 = 0.0_dp
			y1 = muhat_mue(x1,dfdx,lrpar,rpar,lipar,ipar,ierr)
			! check: if muhat < 0, then ground state is pure neutron matter
			if (y1 < 0.0_dp .and. y3 < 0.0_dp) then
				x(i) = 0.0_dp
			else
				x(i) = safe_root(muhat_mue,x1,x3,y1,y3,imax,epsx,epsy,lrpar,rpar,lipar,ipar,ierr)
				if (ierr /= 0) then
					write (*,*) 'unable to converge', rho(i),x1,x3,y1,y3
				end if
			end if
			eps(i) = rho(i)*exp(rpar(3))/clight2
			eps(i) = eps(i)+rho(i)*(Mneutron*(1.0_dp-x(i))+Mproton*x(i))/amu
			P(i) = exp(rpar(4))
			muhat(i) = rpar(5)
			mue(i) = rpar(6)
		end do
	end subroutine find_beta_eq
	
	real(dp) function muhat_mue(x,dfdx,lrpar,rpar,lipar,ipar,ierr)
    ! returns with ierr = 0 if was able to evaluate f and df/dx at x
    ! if df/dx not available, it is okay to set it to 0
    	use constants_def
		use dStar_core_def
		use dStar_eos_def
		use dStar_eos_lib
		
        integer, intent(in) :: lrpar, lipar
        real(dp), intent(in) :: x
        real(dp), intent(out) :: dfdx
        integer, intent(inout), pointer :: ipar(:) ! (lipar)
        real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
        integer, intent(out) :: ierr
		integer :: eos_handle
		real(dp) :: rho, T
		real(dp), dimension(num_dStar_eos_results) :: res
		real(dp) :: muhat, mue
		
		dfdx = 0.0_dp
		eos_handle = ipar(1)
		rho = rpar(1)
		T = rpar(2)
		call eval_core_eos(eos_handle,rho,T,x,res)

		muhat = res(i_mu_n)
		mue = res(i_mu_e)
		muhat_mue = muhat-mue
		
		rpar(3) = res(i_lnE)
		rpar(4) = res(i_lnP)
		rpar(5) = muhat
		rpar(6) = mue

	end function muhat_mue
	
end module beta_equilibrium
