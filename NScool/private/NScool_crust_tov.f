! solves the TOV equations in the crust using pressure as the independent variable

module NScool_crust_tov
    use constants_def

	! the basic variables for the tov solver
	integer, parameter :: tov_radius = 1
	integer, parameter :: tov_baryon = 2
	integer, parameter :: tov_mass = 3
	integer, parameter :: tov_potential = 4
	integer, parameter :: num_tov_variables = 4
    
    integer, parameter :: num_tov_ipar = 0
	
    integer, parameter :: tov_output_step_crust = 1
    integer, parameter :: tov_last_recorded_step = 2
    integer, parameter :: num_tov_rpar = 2
    
	integer, parameter :: tov_default_max_steps = 1000
	real(dp), parameter :: tov_default_max_step_size = 0.0
    real(dp), parameter :: tov_default_starting_step = 1.0e-5
    
contains
    
    subroutine tov_integrate(lgPstart, lgPend, Mcore, Rcore, y, ierr)
        use, intrinsic :: iso_fortran_env, only: error_unit
        real(dp), intent(in) :: lgPstart    ! cgs
        real(dp), intent(in) :: lgPend      ! cgs
        real(dp), intent(in) :: Mcore       ! Msun
        real(dp), intent(in) :: Rcore       ! km
        real(dp), dimension(num_tov_variables), intent(out) :: y
        integer, intent(out) :: ierr
        integer ,dimension(:), allocatable :: iwork
        real(dp), dimension(:), allocatable :: work
        real(dp), dimension(1) :: rtol, atol
        integer, dimension(num_tov_ipar) :: ipar
        integer, dimension(num_tov_rpar) :: rpar
        real(dp) :: lnP, lnPend, h
        integer :: liwork, lwork, itol, lipar, lrpar
        integer :: n, idid, lout, iout
        
		call dop853_work_sizes(num_tov_variables,num_tov_variables,liwork,lwork)
        allocate(iwork(liwork), work(lwork))
        iwork = 0
        iwork(5) = num_tov_variables
        work = 0.0
        
        itol = 0
        iout = 2    ! want dense output
        lout = error_unit
        lipar = num_tov_ipar
        lrpar = num_tov_rpar
        
        rpar(tov_output_step_crust) = 0.1
        rpar(tov_last_recorded_step) = -0.1
        
		y(tov_radius)      = Rcore*1.0e5/length_g
		y(tov_baryon)      = Mcore
		y(tov_mass)        = Mcore
		y(tov_potential)   = 0.0

        n = num_tov_variables
        lnP = lgPstart * ln10
        lnPend = lgPend * ln10
        h = -tov_default_starting_step

		call dop853(n,tov_derivs_crust,lnP,y,lnPend,h,tov_default_max_step_size,tov_default_max_steps, &
			& rtol,atol,itol, tov_solout_crust, iout, work, lwork, iwork, liwork,  &
			&	num_tov_rpar, rpar, num_tov_ipar, ipar, lout, idid)
		deallocate(work, iwork)
    end subroutine tov_integrate
    
    subroutine tov_derivs_crust(n,lnP,y,dy,lrpar,rpar,lipar,ipar,ierr)
    	use dStar_crust_lib

    	integer, intent(in) :: n, lrpar, lipar
    	real(dp), intent(in) :: lnP
    	real(dp), intent(inout) :: y(n)
    	real(dp), intent(out) :: dy(n)
    	real(dp), intent(inout), target :: rpar(lrpar)
    	integer, intent(inout), target :: ipar(lipar)
    	integer, intent(out) :: ierr
    	real(dp) ::  r,r2, r3, a, m, Phi, P, Lnu, Lambda, Hfac, Gfac, lgP, Cv, fourpir2, g
    	real(dp) :: rho, eps, eps_g, rho_g, lgT, nn, np, enu, enu_g, T
        real(dp) :: lgRho, dlgRho, lgNb, dlgNb
        
        ierr = 0
        
        ! everything is in gravitational units
    	P = exp(lnP)
    	r = y(tov_radius)
    	a = y(tov_baryon)
    	m = y(tov_mass)
    	Phi = y(tov_potential)

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

        ! correction factors, see Thorne (1977)
		Lambda = 1.0/sqrt(1.0-2.0*m/r)
		Hfac = 1.0+P/eps_g
		Gfac = (1.0 + fourpi*r3*P/m)
        g = m/r2 * Gfac*Lambda

		dy(tov_radius)      = -P/g/eps_g/Hfac/Lambda
		dy(tov_baryon)      = -P*fourpir2/g/Hfac * rho_g/eps_g
		dy(tov_mass)        = -P*fourpir2/g/Hfac/Lambda
		dy(tov_potential)   = -P/eps_g/Hfac

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
    	real(dp) ::  lnP, lgP, xwant, lgx, r, a, m, phi, p
    	real(dp) :: rho, eps, lgRho, dlgRho, lgNb, dlgNb
	
        ierr = 0
        irtrn = 0
		xwant = rpar(tov_last_recorded_step) - rpar(tov_output_step_crust)
		do while (xwant > x)
			
            lnP = xwant
			r = interp_y(tov_radius, xwant, rwork_y, iwork_y, ierr)
			a = interp_y(tov_baryon, xwant, rwork_y, iwork_y, ierr)
			m = interp_y(tov_mass, xwant, rwork_y, iwork_y, ierr)
			phi = interp_y(tov_potential, xwant, rwork_y, iwork_y, ierr)

			lgP = lnP/ln10 + log10(pressure_g)			! convert P to cgs
            p = exp(lnP)
        	call dStar_crust_get_results(lgP,lgRho,dlgRho,lgNb,dlgNb,ierr)
!             ! our crust model doensn't distinguish between amu and m; small error, will fix soon.
!             eps = 10.0**(lgRho)*clight**2      ! erg cm**-3
!             rho = 10.0**(lgRho)               ! g cm**-3
            
			write (*,'(5(f12.8,tr2),2(es15.8,tr1))') a, m, r*length_g*1.0e-5, 1.0/sqrt(1.0-2.0*m/r), phi, p, 10.0**lgRho

			rpar(tov_last_recorded_step) = rpar(tov_last_recorded_step) - rpar(tov_output_step_crust)
			xwant = rpar(tov_last_recorded_step) - rpar(tov_output_step_crust)
    	end do

    end subroutine tov_solout_crust

end module NScool_crust_tov
