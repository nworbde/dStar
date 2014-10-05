! solves the TOV equations in the crust using pressure as the independent variable

module NScool_crust_tov
    use constants_def

	! the basic variables for the tov solver
	integer, parameter :: tov_radius = 1
	integer, parameter :: tov_baryon = 2
	integer, parameter :: tov_mass = 3
	integer, parameter :: tov_potential = 4
    integer, parameter :: tov_volume = 5
	integer, parameter :: num_tov_variables = 5
    
    integer, parameter :: num_tov_ipar = 0
	
    integer, parameter :: tov_output_step_crust = 1
    integer, parameter :: tov_last_recorded_step = 2
    integer, parameter :: tov_core_mass = 3
    integer, parameter :: tov_core_radius = 4
    integer, parameter :: num_tov_rpar = 4
    
	integer, parameter :: tov_default_max_steps = 1000
	real(dp), parameter :: tov_default_max_step_size = 0.0
    real(dp), parameter :: tov_default_starting_step = 1.0e-5
    
    integer, parameter :: tov_model_expansion_count = 100
    
    type tov_model_type
        ! these are all in gravitational units (c = G = 1; M in units of Msun)
        integer :: nzs
        real(dp) :: core_mass
        real(dp) :: core_radius
        real(dp), allocatable, dimension(:) :: pressure
        real(dp), allocatable, dimension(:) :: radius
        real(dp), allocatable, dimension(:) :: baryon
        real(dp), allocatable, dimension(:) :: mass
        real(dp), allocatable, dimension(:) :: potential
        real(dp), allocatable, dimension(:) :: volume
    end type tov_model_type
    
    type(tov_model_type), target, save :: tov_model
    
contains
    
    subroutine alloc_tov_model(npts, s, ierr)
        integer, intent(in) :: npts
        type (tov_model_type), pointer :: s
        integer, intent(out) :: ierr
        
        allocate(s% pressure(npts), s% radius(npts), s% baryon(npts), s% mass(npts),  &
        &   s% potential(npts), s% volume(npts), stat=ierr)
    end subroutine alloc_tov_model
    
    subroutine free_tov_model(s, ierr)
        type (tov_model_type), pointer :: s
        integer, intent(out) :: ierr
        
        s => tov_model
        deallocate(s% pressure, s% radius, s% baryon, s% mass, s% potential, s% volume, stat=ierr)
    end subroutine free_tov_model
    
    subroutine copy_tov_model(s, snew, ierr)
        type (tov_model_type), pointer :: s, snew
        integer, intent(out) :: ierr
        integer :: n
        ierr = 0

        if (size(snew% pressure) < size(s% pressure)) then
            ierr = -1
            return
        end if
        
        n = s% nzs
        snew% pressure(1:n) = s% pressure(1:n)
        snew% radius(1:n) = s% radius(1:n)
        snew% baryon(1:n) = s% baryon(1:n)
        snew% mass(1:n) = s% mass(1:n)
        snew% potential(1:n) = s% potential(1:n)
        snew% volume(1:n) = s% volume(1:n)
        snew% nzs = s% nzs
        snew% core_mass = s% core_mass
        snew% core_radius = s% core_radius
    end subroutine copy_tov_model
    
    subroutine expand_tov_model(ierr)
        integer, intent(out) :: ierr
        type (tov_model_type), pointer :: s, stmp
        integer :: npts, new_npts
        type(tov_model_type), target :: tmp_model
        
        s => tov_model
        stmp => tmp_model
        
        ierr = 0
        npts = size(s% radius)
        new_npts = npts + tov_model_expansion_count
        call alloc_tov_model(new_npts, stmp, ierr)
        if (ierr /= 0) return
        call copy_tov_model(s, stmp, ierr)
        if (ierr /= 0) return
        call free_tov_model(s, ierr)
        if (ierr /=0) return
        call alloc_tov_model(new_npts, s, ierr)
        if (ierr /= 0) return
        call copy_tov_model(stmp,s,ierr)
        if (ierr /=0) return
        call free_tov_model(stmp, ierr)
    end subroutine expand_tov_model
    
    subroutine tov_integrate(lgPstart, lgPend, Mcore, Rcore, lnP_rez, y, ierr)
        use, intrinsic :: iso_fortran_env, only: error_unit
        use num_lib
        
        real(dp), intent(in) :: lgPstart    ! cgs
        real(dp), intent(in) :: lgPend      ! cgs
        real(dp), intent(in) :: Mcore       ! Msun
        real(dp), intent(in) :: Rcore       ! km
        real(dp), intent(in) :: lnP_rez     ! disired step in lnP
        real(dp), dimension(:), pointer :: y
        integer, intent(out) :: ierr
        integer ,dimension(:), pointer :: iwork => null()
        real(dp), dimension(:), pointer :: work => null()
        real(dp), dimension(1) :: rtol, atol
        integer, dimension(:), pointer :: ipar => null()
        real(dp), dimension(:), pointer :: rpar => null()
        real(dp) :: lnP, lnPend, h
        real(dp) :: phi_correction
        type(tov_model_type), pointer :: s
        integer :: liwork, lwork, itol, lipar, lrpar
        integer :: n, idid, lout, iout, npts
        
        s => tov_model
        
		call dop853_work_sizes(num_tov_variables,num_tov_variables,liwork,lwork)
        allocate(iwork(liwork), work(lwork),ipar(num_tov_ipar), rpar(num_tov_rpar))
        iwork = 0
        iwork(5) = num_tov_variables
        work = 0.0
        
        itol = 0
        rtol = 1.0e-4
        atol = 1.0e-5
        iout = 2    ! want dense output
        lout = error_unit
        lipar = num_tov_ipar
        lrpar = num_tov_rpar
                
		y(tov_radius)      = 0.0
		y(tov_baryon)      = 0.0
		y(tov_mass)        = 0.0
		y(tov_potential)   = 0.5*log(1.0-2.0*Mcore/(Rcore*1.0e5/length_g))
        y(tov_volume)      = 0.0

        n = num_tov_variables
        lnP = lgPstart * ln10 - log(pressure_g)
        lnPend = lgPend * ln10 - log(pressure_g)
        h = -0.1
        rpar(tov_output_step_crust) = lnP_rez
        rpar(tov_last_recorded_step) = lnP + rpar(tov_output_step_crust)
        rpar(tov_core_mass) = Mcore
        rpar(tov_core_radius) = Rcore*1.0e5/length_g

        npts = ceiling(lnP - lnPend)/rpar(tov_output_step_crust) + 1
        
        call alloc_tov_model(npts, s, ierr)
        s% nzs = 0
        s% core_mass = rpar(tov_core_mass)
        s% core_radius = rpar(tov_core_radius)
        
        if (ierr /= 0) return
        
		call dop853(n,tov_derivs_crust,lnP,y,lnPend,h,tov_default_max_step_size,tov_default_max_steps, &
			& rtol,atol,itol, tov_solout_crust, iout, work, lwork, iwork, liwork,  &
			&	num_tov_rpar, rpar, num_tov_ipar, ipar, lout, idid)

		deallocate(work, iwork)
        
        ! now apply the boundary condition to correct phi: exp(phi) = sqrt(1-2m/r)
        phi_correction  = 0.5*log(1.0-2.0*(y(tov_mass)+rpar(tov_core_mass))/(y(tov_radius)+rpar(tov_core_radius))) &
            & -y(tov_potential)
        
        tov_model% potential(1: tov_model% nzs) = tov_model% potential(1: tov_model% nzs) + phi_correction
    end subroutine tov_integrate
    
    subroutine tov_derivs_crust(n,lnP,h,y,dy,lrpar,rpar,lipar,ipar,ierr)
    	use dStar_crust_lib

    	integer, intent(in) :: n, lrpar, lipar
    	real(dp), intent(in) :: lnP, h
    	real(dp), intent(inout) :: y(n)
    	real(dp), intent(out) :: dy(n)
    	real(dp), intent(inout), pointer :: rpar(:)
    	integer, intent(inout), pointer :: ipar(:)
    	integer, intent(out) :: ierr
    	real(dp) ::  r,r2, r3, a, m, Phi, P, Lnu, Lambda, Hfac, Gfac, lgP, Cv, fourpir2, g
    	real(dp) :: rho, eps, eps_g, rho_g, lgT, nn, np, enu, enu_g, T
        real(dp) :: lgRho, dlgRho, lgEps, dlgEps
        
        ierr = 0
        
        ! everything is in gravitational units
    	P = exp(lnP)
    	r = y(tov_radius) + rpar(tov_core_radius)
    	a = y(tov_baryon) + rpar(tov_core_mass)
    	m = y(tov_mass) + rpar(tov_core_mass)
    	Phi = y(tov_potential)

    	lgP = lnP/ln10 + log10(pressure_g)			! convert P to cgs
    	call dStar_crust_get_results(lgP,lgRho,dlgRho,lgEps,dlgEps,ierr)

    	eps = 10.0**(lgEps)	  ! mass-energy density, in g cm**-3
        rho = 10.0**(lgRho)               ! g cm**-3
    	r2 = r*r
    	r3 = r2*r
    	fourpir2 = fourpi*r2

    	! scale to gravitational units
    	rho_g = rho / density_g
    	eps_g = eps / density_g

        ! correction factors, see Thorne (1977)
		Lambda = 1.0/sqrt(1.0-2.0*m/r)
		Hfac = 1.0+P/eps_g
		Gfac = (1.0 + fourpi*r3*P/m)
        g = m/r2 * Gfac*Lambda

		dy(tov_radius)      = -P/g/eps_g/Hfac/Lambda
		dy(tov_baryon)      = -P*fourpir2/g/Hfac * rho_g/eps_g
		dy(tov_mass)        = -P*fourpir2/g/Hfac/Lambda
		dy(tov_potential)   = -P/eps_g/Hfac
        dy(tov_volume)      = fourpir2*Lambda*dy(tov_radius)
        
    end subroutine tov_derivs_crust

    subroutine tov_solout_crust(nr, xold, x, n, y, rwork_y, iwork_y, interp_y, lrpar, rpar, lipar, ipar, irtrn)
        use dStar_crust_lib
        
    	integer, intent(in) :: nr, n, lrpar, lipar
    	real(dp), intent(in) :: xold, x
    	real(dp), intent(inout) :: y(n)
        real(dp), intent(inout), target :: rwork_y(*)
        integer, intent(inout), target :: iwork_y(*)
        integer, intent(inout), pointer :: ipar(:) ! (lipar)
        real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
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
    	real(dp) ::  lnP, lgP, xwant, lgx, r, a, m, phi, p, vol
    	real(dp) :: rho, eps, lgRho, dlgRho, lgEps, dlgEps
	    type(tov_model_type), pointer :: s
        integer :: iz
        
        ierr = 0
        irtrn = 0
		xwant = rpar(tov_last_recorded_step) - rpar(tov_output_step_crust)
		s => tov_model
        
        do while (xwant > x)
			
            iz = s% nzs + 1
            if (iz > size(s% radius)) call expand_tov_model(ierr)
            if (ierr /= 0) exit
            
            lnP = xwant
!             r = interp_y(tov_radius, xwant, rwork_y, iwork_y, ierr) + rpar(tov_core_radius)
!             a = interp_y(tov_baryon, xwant, rwork_y, iwork_y, ierr) + rpar(tov_core_mass)
!             m = interp_y(tov_mass, xwant, rwork_y, iwork_y, ierr) + rpar(tov_core_mass)
!             phi = interp_y(tov_potential, xwant, rwork_y, iwork_y, ierr)
!             vol = interp_y(tov_volume, xwant, rwork_y, iwork_y, ierr)
!
!             lgP = lnP/ln10 + log10(pressure_g)            ! convert P to cgs
!             p = exp(lnP)
!             call dStar_crust_get_results(lgP,lgRho,dlgRho,lgEps,dlgEps,ierr)
            
            s% pressure(iz) = exp(lnP)
            s% radius(iz) = interp_y(tov_radius, xwant, rwork_y, iwork_y, ierr)
            s% baryon(iz) = interp_y(tov_baryon, xwant, rwork_y, iwork_y, ierr)
            s% mass(iz)   = interp_y(tov_mass, xwant, rwork_y, iwork_y, ierr)
            s% potential(iz) = interp_y(tov_potential, xwant, rwork_y, iwork_y, ierr)
            s% volume(iz) = interp_y(tov_volume, xwant, rwork_y, iwork_y, ierr)
            s% nzs = iz
            
!             write (*,'(5(f14.10,tr2),4(es15.8,tr1))') a, m, r*length_g*1.0e-5, 1.0/sqrt(1.0-2.0*m/r), phi,  &
!             &   10.0**lgP, 10.0**lgRho, 10.0**lgEps, vol*length_g**3

			rpar(tov_last_recorded_step) = rpar(tov_last_recorded_step) - rpar(tov_output_step_crust)
			xwant = rpar(tov_last_recorded_step) - rpar(tov_output_step_crust)
            
    	end do
        if (ierr /= 0) irtrn = -1
    end subroutine tov_solout_crust
    
    subroutine tov_write_crust()
        use dStar_crust_lib
    	real(dp) ::  lnP, lgP, xwant, lgx, r, a, m, phi, p, vol
    	real(dp) :: rho, eps, lgRho, dlgRho, lgEps, dlgEps
	    type(tov_model_type), pointer :: s
        integer :: i, ierr

        s => tov_model
        
        do i = 1, s% nzs
            r = s% radius(i) + s% core_radius
            a = s% baryon(i) + s% core_mass
            m = s% mass(i) + s% core_mass
            phi = s% potential(i)
            vol = s% volume(i)

            lgP = log10(s% pressure(i)) + log10(pressure_g)            ! convert P to cgs
            call dStar_crust_get_results(lgP,lgRho,dlgRho,lgEps,dlgEps,ierr)
            if (ierr /= 0) return
            write (*,'(5(f14.10,tr2),4(es15.8,tr1))') a, m, r*length_g*1.0e-5, 1.0/sqrt(1.0-2.0*m/r), phi,  &
            &   10.0**lgP, 10.0**lgRho, 10.0**lgEps, vol*length_g**3
        end do
    end subroutine tov_write_crust

end module NScool_crust_tov
