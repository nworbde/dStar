! solves the TOV equations in the crust using pressure as the independent variable

module NScool_crust_tov
    use constants_def

	! the basic variables for the tov solver
	integer, parameter :: tov_radius = 1
	integer, parameter :: tov_baryon = 2
	integer, parameter :: tov_mass = 3
	integer, parameter :: tov_potential = 4
	integer, parameter :: tov_pressure = 5
	integer, parameter :: num_basic_tov_variables = 5

contains
    
    subroutine tov_derivs_crust(n,lnP,y,dy,lrpar,rpar,lipar,ipar,ierr)
    	use dStar_crust_lib

    	integer, intent(in) :: n, lrpar, lipar
    	real(dp), intent(in) :: lnP
    	real(dp), intent(inout) :: y(n)
    	real(dp), intent(out) :: dy(n)
    	real(dp), intent(inout), target :: rpar(lrpar)
    	integer, intent(inout), target :: ipar(lipar)
    	integer, intent(out) :: ierr
    	real(dp) ::  r,r2, r3, a, m, Phi, P, Lnu, Lambda, Hfac, Gfac, lgP, Cv, fourpir2
    	real(dp) :: rho, eps, eps_g, rho_g, lgT, nn, np, enu, enu_g, T
        
        ierr = 0
        
        ! everything is in gravitational units
    	P = exp(lnP)
    	r = y(tov_radius)
    	a = y(tov_baryon)
    	m = y(tov_mass)
    	Phi = y(tov_potential)
    	P = y(tov_pressure)

    	lgP = lnP/ln10 + log10(pressure_g)			! convert P to cgs
    	call dStar_crust_get_results(lgP,lgRho,dlgRho,lgNb,dlgNb,ierr)
        ! our crust model doensn't distinguish between amu and m; small error, will fix soon.
    	eps = 10.0**(lgRho)*clight**2	  ! erg cm**-3
        rho = 10.0**(lgRho)               ! g cm**-3
    	r2 = r*r
    	r3 = r2*r
    	fourpir2 = fourpi*r2

    	! scale to gravitational units
    	rho_g = rho / density_g
    	eps_g = eps / pressure_g
    	dy = 0.0
		Lambda = 1.0/sqrt(1.0-2.0*m/r)
		Hfac = 1.0+P/eps_g
		Gfac = (1.0 + fourpi*r3*P/m)

		dy(tov_radius) = -P*r2/m/eps_g/Hfac/Lambda**2/Gfac
		dy(tov_baryon) = -P*fourpir2*r2/m/Hfac/Lambda/Gfac * rho_g/eps_g
		dy(tov_mass) = -fourpir2*r2*P/m/Hfac/Lambda**2/Gfac
		dy(tov_potential) = -P/eps_g/Hfac
		dy(tov_pressure) = P

    end subroutine tov_derivs_crust

    subroutine tov_solout_crust(nr, xold, x, n, y, rwork_y, iwork_y, interp_y, lrpar, rpar, lipar, ipar, irtrn)
        use dStar_crust_lib
        
    	integer, intent(in) :: nr, n, lrpar, lipar
    	real(dp), intent(in) :: xold, x
    	real(dp), intent(inout) :: y(n)
    	real(dp), intent(inout), target :: rpar(lrpar), rwork_y(*)
    	integer, intent(inout), target :: ipar(lipar), iwork_y(*)
    	integer, intent(out) :: irtrn
    	interface
    		double precision function interp_y(i, s, rwork_y, iwork_y, ierr)
            integer, intent(in) :: i ! result is interpolated approximation of y(i) at x=s.
            double precision, intent(in) :: s ! interpolation x value (between xold and x).
            double precision, intent(inout), target :: rwork_y(*)
            integer, intent(inout), target :: iwork_y(*)
            integer, intent(out) :: ierr
         end function interp_y
    	end interface ! 
    	integer :: ierr, i
    	real(dp) ::  lgP, xwant, lgx, r, a, m, phi, p
    	real(dp) :: rho, eps
	
        ierr = 0
        irtrn = 0
		xwant = rpar(tov_last_recorded_step) - rpar(tov_output_step_crust)
		do while (xwant > x)
			
			r = interp_y(tov_radius, xwant, rwork_y, iwork_y, ierr)
			a = interp_y(tov_baryon, xwant, rwork_y, iwork_y, ierr)
			m = interp_y(tov_mass, xwant, rwork_y, iwork_y, ierr)
			phi = interp_y(tov_potential, xwant, rwork_y, iwork_y, ierr)
			p = interp_y(tov_pressure, xwant, rwork_y, iwork_y, ierr)

			lgP = log10(y(tov_pressure))+log10(pressure_g)			! convert P to cgs
        	call dStar_crust_get_results(lgP,lgRho,dlgRho,lgNb,dlgNb,ierr)
!             ! our crust model doensn't distinguish between amu and m; small error, will fix soon.
!             eps = 10.0**(lgRho)*clight**2      ! erg cm**-3
!             rho = 10.0**(lgRho)               ! g cm**-3
            
			write (*,'(4(es15.8,tr1))') a, m, r*length_g, phi, p*pressure_g, 10.0**lgRho

			rpar(tov_last_recorded_step) = rpar(tov_last_recorded_step) - rpar(tov_output_step_crust)
			xwant = rpar(tov_last_recorded_step) - rpar(tov_output_step_crust)
    	end do

    end subroutine tov_solout_crust


end module NScool_crust_tov
