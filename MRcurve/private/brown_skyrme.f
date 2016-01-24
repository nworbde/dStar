module brown_skyrme
	use constants_def
	integer, parameter :: skyrme_matter = 1, skyrme_neutron = 2, skyrme_sym = 3

	type skyrme
		integer :: model
		real(dp) :: gamma
		real(dp) :: a, b, c, d, e
	end type skyrme

	type(skyrme), target :: skyrme_eos(3)
	
contains
	
	subroutine initialize_brown_skyrme(ierr)
		integer, intent(out) :: ierr
		type(skyrme), pointer :: s
		s => skyrme_eos(skyrme_matter)
		s% model = skyrme_matter
		s% gamma = 1.25000
		s% a = -807.32
		s% b = 914.73
		s% c = 75.01
		s% d = -58.95
		s% e = 312.50

		s => skyrme_eos(skyrme_neutron)
		s% model = skyrme_sym
		s% gamma = 1.25000
		s% a = -377.46
		s% b = 406.07
		s% c = 118.98
		s% d = -5.37
		s% e = 138.73

		s => skyrme_eos(skyrme_sym)
		s% model = skyrme_sym
		s% gamma = 1.25000
		s% a = 429.86
		s% b = -508.66
		s% c = 43.98
		s% d = 53.59 
		s% e = -173.77		
		
		ierr = 0
	end subroutine initialize_brown_skyrme
	
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
			x1 = 0.0_dp
			y1 = muhat_mue(x1,dfdx,lrpar,rpar,lipar,ipar,ierr)
			x3 = 0.5_dp
			y3 = muhat_mue(x1,dfdx,lrpar,rpar,lipar,ipar,ierr)
			x(i) = safe_root(muhat_mue,x1,x3,y1,y3,imax,epsx,epsy,lrpar,rpar,lipar,ipar,ierr)
			if (ierr /= 0) then
				write (*,*) 'unable to converge', rho(i),x1,x3,y1,y3
				cycle
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
		muhat_mue = max(muhat-mue,0.0_dp)
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
		eps = rho*(Eint+ x*Mp_n + (1.0-x)*Mn_n)
		P = Pskyrme(matter,rho)+xfac*Pskyrme(sym,rho)
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
		
	end subroutine get_skyrme_eos
	
end module brown_skyrme
