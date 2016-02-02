module brown_skyrme
	use constants_def
	use dStar_core_def
contains
	
	subroutine find_beta_eq(rho,x,eps,P,muhat,mue,cs2,ierr)
		use constants_def
		use dStar_core_def
		use num_lib
		
		real(dp), dimension(:), intent(in) :: rho
		real(dp), dimension(:), intent(out) :: x,eps,P,muhat,mue,cs2
		integer, intent(out) :: ierr
		
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
		
		lipar = 1
		allocate(ipar(lipar))
		ipar(1) = 6
		
		lrpar = ipar(1)
		allocate(rpar(lrpar))
		
		do i = 1, Ntab
			rpar(1) = rho(i)
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
			eps(i) = rpar(2)
			P(i) = rpar(3)
			muhat(i) = rpar(4)
			mue(i) = rpar(5)
			cs2(i) = rpar(6)
		end do
	end subroutine find_beta_eq
	
	real(dp) function muhat_mue(x,dfdx,lrpar,rpar,lipar,ipar,ierr)
    ! returns with ierr = 0 if was able to evaluate f and df/dx at x
    ! if df/dx not available, it is okay to set it to 0
    	use constants_def
		use dStar_core_def
		
        integer, intent(in) :: lrpar, lipar
        real(dp), intent(in) :: x
        real(dp), intent(out) :: dfdx
        integer, intent(inout), pointer :: ipar(:) ! (lipar)
        real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
        integer, intent(out) :: ierr
		real(dp) :: rho, eps, P, muhat, mue, cs2
		
		dfdx = 0.0_dp
		rho = rpar(1)
		call get_skyrme_eos(rho,x,eps,P,muhat,mue,cs2,ierr)
		if (ierr /= 0) then
			write (*,*) 'unable to compute EOS'
			muhat_mue = 0.0_dp
			return
		end if
		muhat_mue = muhat-mue
		rpar(2) = eps
		rpar(3) = P
		rpar(4) = muhat
		rpar(5) = mue
		rpar(6) = cs2
	end function muhat_mue
	
	subroutine get_skyrme_eos(rho,x,eps,P,muhat,mue,cs2,ierr)
		! rho := nucleon density (fm**-3)
		! x := proton fraction
		! eps := energy density (MeV/fm**3)
		! P := pressure (MeV/fm**3)
		! cs2 := square sound speed (c=1)
		real(dp), intent(in) :: rho,x
		real(dp), intent(out) :: eps,P,muhat,mue,cs2
		integer, intent(out) :: ierr
		type(skyrme), pointer :: matter, sym
		real(dp) :: xfac, dEdlnrho, Eint, Mhat
		
		Mhat = Mn_n - Mp_n
		
		ierr = 0
		matter => skyrme_eos(skyrme_matter)
		sym => skyrme_eos(skyrme_sym)
		xfac = (1.0-2.0*x)**2
		Eint = Eskyrme(matter,rho) + xfac*Eskyrme(sym,rho)
		eps = rho*(Eint+ x*Mp_n + (1.0-x)*Mn_n) + eps_e(rho,x)
		P = Pskyrme(matter,rho)+xfac*Pskyrme(sym,rho) + P_e(rho,x)
		muhat = 4.0*(1.0-2.0*x)*Eskyrme(sym,rho)+Mhat
		mue = mu_e(rho,x)
		dEdlnrho = rhodEdrho(matter,rho)+xfac*rhodEdrho(sym,rho)
		cs2 = rho*dEdlnrho/(Eint + dEdlnrho)
		
	contains
		function Eskyrme(s,rho)
			type(skyrme), pointer :: s
			real(dp), intent(in) :: rho
			real(dp) :: Eskyrme
	        Eskyrme = s% a*rho + s% b*rho**s% gamma + s% c*rho**twothird +  &
	        	& s% d*rho**fivethird + s% e*rho**3
		end function Eskyrme
		
	    function rhodEdrho(s,rho)
			type(skyrme), pointer :: s
			real(dp), intent(in) :: rho
			real(dp) :: rhodEdrho
	        rhodEdrho =  s% a*rho + s% b*s% gamma*rho**s% gamma +  &
	        	& twothird*s% c*rho**twothird + fivethird*s% d*rho**fivethird &
	        	& +3.0*s% e*rho**3
		end function rhodEdrho
	    
		function Pskyrme(s,rho)
	        type(skyrme), pointer :: s
			real(dp), intent(in) :: rho
			real(dp) :: Pskyrme
			Pskyrme = rho*(s% a*rho + s% b*s% gamma*rho**s% gamma +  &
	        	& twothird*s% c*rho**twothird + fivethird*s% d*rho**fivethird &
	        	& +3.0*s% e*rho**3)
		end function Pskyrme
    
	    function dPdrho(s,rho)
			type(skyrme), pointer :: s
			real(dp), intent(in) :: rho
			real(dp) :: dPdrho
	        dPdrho = 2*s% a*rho + s% b*s% gamma*(s% gamma+1.0)*rho**s% gamma  &
	        	& + s% c*twothird*fivethird*rho**twothird  &
	        	& + s% d*fivethird*seventhird*rho**fivethird + 12.0*s% e*rho**3
		end function dPdrho
		
		function mu_e(rho,x)
			real(dp), intent(in) :: rho,x
			real(dp) :: mu_e
			mu_e = (threepisquare*x*rho)**onethird*hbarc_n
		end function mu_e
		
		function P_e(rho,x)
			real(dp), intent(in) :: rho,x
			real(dp) :: P_e
			P_e = 0.25*(x*rho)*(threepisquare*x*rho)**onethird*hbarc_n
		end function P_e
		
		function eps_e(rho,x)
			real(dp), intent(in) :: rho,x
			real(dp) :: eps_e
			eps_e = 0.75*(x*rho)*(threepisquare*x*rho)**onethird*hbarc_n
		end function eps_e
		
	end subroutine get_skyrme_eos
	
end module brown_skyrme
