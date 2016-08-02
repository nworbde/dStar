! solves the TOV equations using pressure as the independent variable

module dStar_core_tov
    use constants_def

	! the basic variables for the tov solver
	integer, parameter :: core_radius = 1
	integer, parameter :: core_baryon = 2
	integer, parameter :: core_mass = 3
	integer, parameter :: core_potential = 4
	integer, parameter :: num_core_variables = 4
    
    integer, parameter :: num_core_ipar = 0
	
    integer, parameter :: core_output_step_crust = 1
    integer, parameter :: core_last_recorded_step = 2
    integer, parameter :: num_core_rpar = 2
    
	integer, parameter :: core_default_max_steps = 1000
	real(dp), parameter :: core_default_max_step_size = 0.0
    real(dp), parameter :: core_default_starting_step = 1.0e-5
    
    integer, parameter :: core_model_expansion_count = 100
    
	real(dp), parameter :: lgP_blend_max = 0.35	! MeV fm**-3
	real(dp), parameter :: lgP_blend_min = -0.80	! MeV fm**-3
	
    type core_structure_type
        ! these are all in gravitational units (c = G = 1; M in units of Msun)
        integer :: nzs
        real(dp), allocatable, dimension(:) :: pressure
        real(dp), allocatable, dimension(:) :: radius
        real(dp), allocatable, dimension(:) :: baryon
        real(dp), allocatable, dimension(:) :: mass
        real(dp), allocatable, dimension(:) :: potential
    end type core_structure_type
    
    type(core_structure_type), target, save :: core_structure
    
contains
    
    subroutine alloc_core_model(npts, s, ierr)
        integer, intent(in) :: npts
        type (core_structure_type), pointer :: s
        integer, intent(out) :: ierr
        
        allocate(s% pressure(npts), s% radius(npts), s% baryon(npts), s% mass(npts),  &
        &   s% potential(npts), stat=ierr)
    end subroutine alloc_core_model
    
    subroutine free_core_model(s, ierr)
        type (core_structure_type), pointer :: s
        integer, intent(out) :: ierr
        
        s => core_structure
        deallocate(s% pressure, s% radius, s% baryon, s% mass, s% potential, stat=ierr)
    end subroutine free_core_model
    
    subroutine copy_core_model(s, snew, ierr)
        type (core_structure_type), pointer :: s, snew
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
        snew% nzs = s% nzs
    end subroutine copy_core_model
    
    subroutine expand_core_model(ierr)
        integer, intent(out) :: ierr
        type (core_structure_type), pointer :: s, stmp
        integer :: npts, new_npts
        type(core_structure_type), target :: tmp_model
        
        s => core_structure
        stmp => tmp_model
        
        ierr = 0
        npts = size(s% radius)
        new_npts = npts + core_model_expansion_count
        call alloc_core_model(new_npts, stmp, ierr)
        if (ierr /= 0) return
        call copy_core_model(s, stmp, ierr)
        if (ierr /= 0) return
        call free_core_model(s, ierr)
        if (ierr /=0) return
        call alloc_core_model(new_npts, s, ierr)
        if (ierr /= 0) return
        call copy_core_model(stmp,s,ierr)
        if (ierr /=0) return
        call free_core_model(stmp, ierr)
    end subroutine expand_core_model
    
    subroutine core_integrate(lgPstart, lgPend, lnP_rez, y, ierr)
        use, intrinsic :: iso_fortran_env, only: error_unit
        use num_lib
        
        real(dp), intent(in) :: lgPstart    ! MeV fm**-3
        real(dp), intent(in) :: lgPend      ! cgs
        real(dp), intent(in) :: lnP_rez     ! desired step in lnP
        real(dp), dimension(:), pointer :: y	! solution vector
        integer, intent(out) :: ierr
        integer ,dimension(:), pointer :: iwork => null()
        real(dp), dimension(:), pointer :: work => null()
        real(dp), dimension(1) :: rtol, atol
        integer, dimension(:), pointer :: ipar => null()
        real(dp), dimension(:), pointer :: rpar => null()
        real(dp) :: lnP, lnPend, h
        real(dp) :: phi_correction
		real(dp) :: lgRho,lgEps,eps,rho,rho_g,eps_g
        type(core_structure_type), pointer :: s
        integer :: liwork, lwork, itol, lipar, lrpar
        integer :: n, idid, lout, iout, npts
        
        s => core_structure
        
		call dop853_work_sizes(num_core_variables,num_core_variables,liwork,lwork)
        allocate(iwork(liwork), work(lwork),ipar(num_core_ipar), rpar(num_core_rpar))
        iwork = 0
        iwork(5) = num_core_variables
        work = 0.0
        
        itol = 0
        rtol = 1.0e-5
        atol = 1.0e-5
        iout = 2    ! want dense output.. eventually turn off?
        lout = error_unit
        lipar = num_core_ipar
        lrpar = num_core_rpar
                
		! initial values: start with small 10 m sphere of constant density
		call core_get_EOS(lgPstart,lgRho,lgEps,ierr)
		if (ierr /= 0) return
		
    	eps = 10.0**(lgEps)	  ! mass-energy density, in MeV fm**-3
        rho = 10.0**(lgRho)   ! fm**-3
		
		write(error_unit,'(a,2(e11.4))') 'integrating with central densities',rho*amu*density_n,eps*mass_n*density_n

    	! scale to gravitational units
    	rho_g = rho * amu*density_n/density_g
    	eps_g = eps * mass_n*density_n/density_g
		
		y(core_radius)      = 3.0e3_dp/length_g
		y(core_baryon)      = onethird*fourpi*rho_g*y(core_radius)**3
		y(core_mass)        = y(core_baryon)*eps_g/rho_g
		y(core_potential)   = 0.0

        n = num_core_variables
        lnP = lgPstart * ln10 + log(pressure_n/pressure_g)
        lnPend = lgPend * ln10 + log(pressure_n/pressure_g)
		! make the initial stepsize very small because dP/dr is small at center.
        h = -1.0e-4
        rpar(core_output_step_crust) = lnP_rez
        rpar(core_last_recorded_step) = lnP + rpar(core_output_step_crust)

        npts = ceiling(lnP - lnPend)/rpar(core_output_step_crust) + 1
        
		if (allocated(s% pressure)) then
			call free_core_model(s, ierr)
		end if
        call alloc_core_model(npts, s, ierr)
        s% nzs = 0
        
        if (ierr /= 0) return
        
		call dop853(n,core_derivs_crust,lnP,y,lnPend,h,core_default_max_step_size, &
			& core_default_max_steps, &
			& rtol,atol,itol, core_solout_crust, iout, work, lwork, iwork, liwork,  &
			& num_core_rpar, rpar, num_core_ipar, ipar, lout, idid)

		deallocate(work, iwork)
        
        ! now apply the boundary condition to correct phi: exp(phi) = sqrt(1-2m/r)
        phi_correction  = 0.5*log(1.0-2.0*(y(core_mass))/(y(core_radius))) &
            & -y(core_potential)
        
        core_structure% potential(1: core_structure% nzs) =  &
        	& core_structure% potential(1: core_structure% nzs) + phi_correction
    end subroutine core_integrate
    
	subroutine core_get_EOS(lgP,lgRho,lgEps,ierr)
    	use dStar_crust_lib
		use dStar_core_lib, only : dStar_core_get_results
		use dStar_core_mod
		real(dp), intent(in) :: lgP	! MeV fm**-3
		real(dp), intent(out) :: lgRho, lgEps	! fm**-3, MeV fm**-3
		integer, intent(out) :: ierr
        real(dp) :: dlgRho, dlgEps
		real(dp) :: lgRhon, dlgRhon, lgEpsn, dlgEpsn
		real(dp) :: lgPc, lgRhoc, dlgRhoc, lgEpsc, dlgEpsc, Xprot, dXprot
		real(dp) :: Zfac, fac
		
		! get the pressure in cgs for the crust EOS
		lgPc = lgP + log10(pressure_n)
		
		! now determine our regime and get the appropriate coefficients
		if (lgP <= lgP_blend_max) then
			! we shall need the crust coefficients
			call dStar_crust_get_eos(lgPc,lgRhoc,dlgRhoc,lgEpsc,dlgEpsc,ierr)
			! convert to nuclear units
			lgRhoc = lgRhoc - log10(amu*density_n)
			lgEpsc = lgEpsc - log10(mass_n*density_n)
		end if
		if (lgP > lgP_blend_min) then
			call dStar_core_get_results(lgP,lgRhon,dlgRhon,lgEpsn,dlgEpsn,Xprot,dXprot,ierr)
		end if
		
		! now determine regime and blend if necessary
		if (lgP <= lgP_blend_min) then
			lgRho = lgRhoc
			lgEps = lgEpsc
		else if (lgP <= lgP_blend_max .and. lgP > lgP_blend_min) then
			Zfac = (lgP-lgP_blend_min)/(lgP_blend_max-lgP_blend_min)
			fac = 0.5_dp*(1.0_dp-cos(pi*Zfac))
			lgRho = fac*lgRhon + (1.0-fac)*lgRhoc
			lgEps = fac*lgEpsn + (1.0-fac)*lgEpsc
		else
			lgRho = lgRhon
			lgEps = lgEpsn
		end if
	end subroutine core_get_EOS
	
    subroutine core_derivs_crust(n,lnP,h,y,dy,lrpar,rpar,lipar,ipar,ierr)
    	integer, intent(in) :: n, lrpar, lipar
    	real(dp), intent(in) :: lnP, h
    	real(dp), intent(inout) :: y(:)
    	real(dp), intent(out) :: dy(:)
    	real(dp), intent(inout), pointer :: rpar(:)
    	integer, intent(inout), pointer :: ipar(:)
    	integer, intent(out) :: ierr
    	real(dp) :: r,r2, r3, a, m, Phi, P
    	real(dp) :: Lambda, Hfac, Gfac, lgP, Cv, fourpir2, g
    	real(dp) :: lgRho, lgEps, eps_g, rho_g
        
        ierr = 0
        
        ! everything is in gravitational units
    	P = exp(lnP)
    	r = y(core_radius)
    	a = y(core_baryon)
    	m = y(core_mass)
    	Phi = y(core_potential)

		! convert P to MeV fm**-3
    	lgP = lnP/ln10 + log10(pressure_g/pressure_n)
		call core_get_EOS(lgP,lgRho,lgEps,ierr)
		
		! convert to gravitational units
		rho_g = 10.0**lgRho * amu*density_n / density_g
		eps_g = 10.0**lgEps * mass_n*density_n / density_g
		
    	r2 = r*r
    	r3 = r2*r
    	fourpir2 = fourpi*r2

        ! correction factors, see Thorne (1977)
		Lambda = 1.0/sqrt(1.0-2.0*m/r)
		Hfac = 1.0+P/eps_g
		Gfac = (1.0 + fourpi*r3*P/m)
        g = m/r2 * Gfac*Lambda

		dy(core_radius)      = -P/g/eps_g/Hfac/Lambda
		dy(core_baryon)      = -P*fourpir2/g/Hfac * rho_g/eps_g
		dy(core_mass)        = -P*fourpir2/g/Hfac/Lambda
		dy(core_potential)   = -P/eps_g/Hfac
        
    end subroutine core_derivs_crust

    subroutine core_solout_crust(nr, xold, x, n, y, rwork_y, iwork_y, interp_y, lrpar, rpar, lipar, ipar, irtrn)
    	integer, intent(in) :: nr, n, lrpar, lipar
    	real(dp), intent(in) :: xold, x
    	real(dp), intent(inout) :: y(:)
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
	    type(core_structure_type), pointer :: s
        integer :: iz
        
        ierr = 0
        irtrn = 0
		xwant = rpar(core_last_recorded_step) - rpar(core_output_step_crust)
		s => core_structure
        
        do while (xwant > x)
			
            iz = s% nzs + 1
            if (iz > size(s% radius)) call expand_core_model(ierr)
            if (ierr /= 0) exit
            
            lnP = xwant
            
            s% pressure(iz) = exp(lnP)
            s% radius(iz) = interp_y(core_radius, xwant, rwork_y, iwork_y, ierr)
            s% baryon(iz) = interp_y(core_baryon, xwant, rwork_y, iwork_y, ierr)
            s% mass(iz)   = interp_y(core_mass, xwant, rwork_y, iwork_y, ierr)
            s% potential(iz) = interp_y(core_potential, xwant, rwork_y, iwork_y, ierr)
            s% nzs = iz
            
!             write (*,'(5(f14.10,tr2),4(es15.8,tr1))') a, m, r*length_g*1.0e-5, 1.0/sqrt(1.0-2.0*m/r), phi,  &
!             &   10.0**lgP, 10.0**lgRho, 10.0**lgEps, vol*length_g**3

			rpar(core_last_recorded_step) = rpar(core_last_recorded_step) - rpar(core_output_step_crust)
			xwant = rpar(core_last_recorded_step) - rpar(core_output_step_crust)
            
    	end do
        if (ierr /= 0) irtrn = -1
    end subroutine core_solout_crust
    
    subroutine core_write_crust(prefix,eos)
		character(len=*), intent(in) :: prefix,eos
    	real(dp) ::  lnP, lgP, xwant, lgx, r, a, m, phi, p, vol
    	real(dp) :: rho, eps, lgRho, lgEps
	    type(core_structure_type), pointer :: s
        integer :: i, ierr, unitno
		integer :: M_id
		character(len=256) :: outfile
		
        s => core_structure
        
		M_id = int(s% mass(s% nzs)*1000)
		write(outfile,'(a,"/",a,a,"_",i4.4)') trim(prefix),'profile_',trim(eos),M_id
	
		open(newunit=unitno,file=trim(outfile),status='unknown',action='write')
		write (unitno,'(5(a14,tr2),3(a15,tr1))') 'baryon (Msun)','mass (Msun)', &
		& 'radius (km)','Lambda','Phi','pressure','rho','epsilon'
		
        do i = 1, s% nzs
            r = s% radius(i)
            a = s% baryon(i)
            m = s% mass(i)
            phi = s% potential(i)

			! convert P to MeV fm**-3
	    	lgP = log10(s% pressure(i)) + log10(pressure_g/pressure_n)
			call core_get_EOS(lgP,lgRho,lgEps,ierr)
			
            if (ierr /= 0) return
            write (unitno,'(5(f14.10,tr2),3(es15.8,tr1))') a, m, r*length_g*1.0e-5, 1.0/sqrt(1.0-2.0*m/r), phi,  &
            &   10.0**lgP * pressure_n, 10.0**lgRho * amu*density_n, 10.0**lgEps * mass_n*density_n
        end do
		close(unitno)
    end subroutine core_write_crust

end module dStar_core_tov
