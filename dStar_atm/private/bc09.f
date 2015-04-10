module bc09
    
    ! pointers for the photosphere routine
    integer, parameter :: igrav = 1
    integer, parameter :: itau = 2
    integer, parameter :: iTeff = 3
    integer, parameter :: iPph = 4
    integer, parameter :: iKph = 5
    integer, parameter :: number_photosphere_rpar = 5
    integer, parameter :: ihandle = 1
    integer, parameter :: ichem_id = 2
    integer, parameter :: number_photosphere_ipar = 2

contains
    
	subroutine do_get_bc09_Teff(grav, Plight, Tb, Teff, flux)
		use constants_def
		real(dp), intent(in) :: grav	! surface gravity, in the local frame
		real(dp), intent(in) :: Plight	! pressure at which layer of light elements terminates
		real(dp), intent(in), dimension(:) :: Tb	! temperature at a base column
		real(dp), intent(out), dimension(:) :: Teff, flux	! effective temperature and flux
		real(dp) :: eta, g14
		real(dp), dimension(size(Tb)) :: Tb9, Teff6_4
        
        ! make a very dense table of Tb(Teff); then interpolate to get Teff(Tb)
        integer, parameter ::  size_tab = 4*size(Tb)
        real(dp), dimension(:), allocatable :: tabTb9, tabTeff_4
        
        allocate(tabTb9(size_tab),tabTeff_4(size_tab))
        
        ! compute dense table
        
        do i = 1, size_tab
		tabTeff_4(i) = 1.
		! get Pph(Teff)
		call find_photospheric_pressure(tabTeff_4(i),grav,tau,Pphoto,eos_handle,ierr) 
		write(*,*) tabTeff_4(i), Pphoto           
        end do
        
        ! interpolate from dense table to get finished product
        
        deallocate(tabTb9,tabTeff_4)
    end subroutine do_get_bc09_Teff
    
    subroutine find_photospheric_pressure(Teff,grav,tau,Pphoto,eos_handle,ierr)
        use constants_def
        use nucchem_def, only : nuclide_not_found
    	use nucchem_lib, only : get_nuclide_index
    	use num_lib
        
        real(dp), intent(in) :: Teff,grav,tau
        real(dp), intent(out) :: Pphoto
        integer, intent(in) :: eos_handle
        integer, intent(out) :: ierr
        real(dp) :: root_ph,rho,dfdrho
        integer, pointer :: ipar(:) => null() ! (lipar)
        real(dp), pointer :: rpar(:) => null()  ! (lrpar)
        real(dp) :: rhoph_guess, kap_th
        real(dp) :: x1, x3, y1, y3, epsx, epsy
        integer :: imax
        
        ierr = 0
        allocate(ipar(number_photosphere_ipar), rpar(number_photosphere_rpar))

        ! hardwire photosphere compostion to pure He
        ipar(ichem_id) = get_nuclide_index('he4')
        if (ipar(ichem_id) == nuclide_not_found) then
            ierr = nuclide_not_found
            return
        end if
                
        rpar(igrav) = grav
        rpar(iTeff) = Teff
        rpar(itau) = tau
        ipar(ihandle) = eos_handle

     	! parameters for root find
        imax = 20
        epsx = 1.0d-8
        epsy = 1.0d-8
        
        ! use initial guess with ideal gas pressure and thomson scattering
        kap_th = 8.0_dp*onethird*pi*(electroncharge**2/Melectron/clight2)**2 * (0.5)/amu        
        rhoph_guess = 2.0*onethird*gravity/kap_th/(boltzmann*Teff)*amu
        ! brackets for root find 
        x1 = 1.0d2 	! lower bound [g/cm^3]
        x3 = 1.0d9 	! upper bound [g/cm^3]
        y1 = photosphere(x1,dfdrho,lrpar,rpar,lipar,ipar,ierr)	! y1=Pph(x1)
        y1 = rpar(iPph)
        y3 = photosphere(x3,dfdrho,lrpar,rpar,lipar,ipar,ierr)	! y3=Pph(x3)
        y3 = rpar(iPph)
        if (ierr /= 0) then
        	write(*,*) 'unable to create brackets for photosphere root find'
        	return
        end if

		root_ph = safe_root_with_initial_guess(photosphere,rhoph_guess,x1,x3,y1,y3 &
            &   imax,epsx,epsy,lrpar,rpar,lipar,ipar,ierr)
            if (ierr /= 0) then
                write(*,*) 'unable to converge', rhoph_guess, x1, x3, y1, y3
                cycle
            end if
        Pphoto = rpar(iPph)
        
        deallocate(ipar, rpar)
    end subroutine find_photospheric_pressure    

    real(dp) function photosphere(rho, dfdrho, lrpar, rpar, lipar, ipar, ierr)
       ! returns with ierr = 0 if was able to evaluate f and df/dx at x
       ! if df/dx not available, it is okay to set it to 0
       use constants_def
       use nucchem_def
       use nucchem_lib
       use dStar_eos_lib
       use conductivity_lib
       
       integer, intent(in) :: lrpar, lipar
       real(dp), intent(in) :: rho
       real(dp), intent(out) :: dfdrho
       integer, intent(inout), pointer :: ipar(:) ! (lipar)
       real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
       integer, intent(out) :: ierr
       real(dp) :: gravity, tau_ph, Teff, P, kappa
       integer :: eos_handle
       type(composition_info_type) :: ionic
       type(conductivity_components) :: K
       real(dp), dimension(num_dStar_eos_results) :: res
       integer :: phase, ncharged
       integer, dimension(1) :: charged_ids, chem_ids
	   real(dp),dimension(1) :: Y,Yion
	   
       ierr = 0
       
       ! check inputs
       if (rho < 0.0) then
           ierr = -9
           return
       end if
       
       gravity = rpar(igrav)
       tau_ph = rpar(itau)
       Teff = rpar(iTeff)
       eos_handle = ipar(ihandle)
	   chem_ids = ipar(ichem_id)

       ! single species only
       Y(1) = 1.0/nuclib% A(chem_ids(1))

       call compute_composition_moments(1,chem_ids,Y,ionic,Xsum,ncharged, charged_ids, Yion,  &
   			& exclude_neutrons = .TRUE.)
       call eval_crust_eos(eos_handle,rho,Teff,ion,ncharged,charged_ids,Yion, &
       		&   res,phase,use_default_nuclear_size)
       
       P = exp(res(i_lnP))
       rpar(iPph) = P
       get_thermal_conductivity(rho,Teff,chi,Gamma,eta,ionic,K,cond_use_only_kap)
       kappa = 4.0*onethird*arad*clight*Teff**3/rho/K(icond_kap)
       rpar(iKph) = kappa
	   dfdrho = 0.0
       
       return P - 2.0*onethird*gravity/kappa
    end function photosphere

end module bc09
