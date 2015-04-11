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
    
    ! error codes
    integer, parameter :: negative_photosphere_gas_pressure = -2
    integer, parameter :: negative_photosphere_density = -3

contains
    
	subroutine do_get_bc09_Teff(grav, Plight, Tb, Teff, flux)
		use constants_def
		real(dp), intent(in) :: grav	! surface gravity, in the local frame
		real(dp), intent(in) :: Plight	! pressure at which layer of light elements terminates
		real(dp), intent(in), dimension(:) :: Tb	! temperature at a base column
		real(dp), intent(out), dimension(:) :: Teff, flux	! effective temperature and flux
		real(dp) :: eta, g14
        integer ::  size_tab ! = 4*size(Tb)
        real(dp), dimension(:), allocatable :: tabTb9, tabTeff, tabTeff6_4
        integer :: i
        
        ! make a very dense table of Tb(Teff); then interpolate to get Teff(Tb)        
        size_tab = 4*size(Tb)
        allocate(tabTb9(size_tab),tabTeff(size_tab),tabTeff6_4(size_tab))
!         tau = ?
!         Teff = ?
        ! compute dense table
        do i = 1, size_tab
            ! get Pph(Teff)
!            call find_photospheric_pressure(Teff,grav,tau,Pphoto,eos_handle,ierr) 
!		write(*,*) tabTeff_4(i), Pphoto           
        end do
        
        ! interpolate from dense table to get finished product
        
        deallocate(tabTb9,tabTeff6_4)
    end subroutine do_get_bc09_Teff
    
    subroutine find_photospheric_pressure(Teff,grav,tau,Pphoto,eos_handle,ierr)
        use constants_def
        use nucchem_def, only : nuclide_not_found, nuclib
    	use nucchem_lib, only : get_nuclide_index
    	use num_lib
        integer, parameter :: default_maximum_iterations_photosphere = 20
        real(dp), parameter :: default_tolerance_photosphere_rho = 1.0e-8_dp
        real(dp), parameter :: default_tolerance_photosphere_condition = 1.0e-8_dp
        real(dp), intent(in) :: Teff,grav,tau
        real(dp), intent(out) :: Pphoto
        integer, intent(in) :: eos_handle
        integer, intent(out) :: ierr
        real(dp) :: root_ph,rho,dfdrho
        integer, pointer :: ipar(:) => null() ! (lipar)
        real(dp), pointer :: rpar(:) => null()  ! (lrpar)
        integer :: i, lrpar, lipar
        real(dp) :: sigma_Th
        real(dp) :: fallback_Pphoto
        real(dp) :: rho_guess, kappa_Th, A, Z, Pgas
        real(dp) :: rho1, rho3, drho, y1, y3
        real(dp) :: eps_rho, eps_ph
        integer :: maximum_iterations
        
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
        lrpar = number_photosphere_rpar
        lipar = number_photosphere_ipar
        
        ! use initial guess with ideal gas pressure and thomson scattering
        A = nuclib% A(ipar(ichem_id))
        Z = nuclib% Z(ipar(ichem_id))
        sigma_Th = 8.0_dp*onethird*pi*(electroncharge**2/Melectron/clight2)**2
        kappa_Th = sigma_Th* Z/A
        fallback_Pphoto = 2.0_dp*onethird*arad*Teff**4   
        Pgas = twothird*grav/kappa_th - onethird*arad*Teff**4
        if (Pgas < 0.0_dp) then
            ierr = negative_photosphere_gas_pressure
            Pphoto = fallback_Pphoto
            return
        end if
        rho_guess = Pgas*amu*A/(Z+1.0)/(boltzmann*Teff)
        ! get brackets for root find
        drho = 0.1_dp
        maximum_iterations = 10
        call look_for_brackets(rho_guess, drho, rho1, rho3, photosphere, y1, y3, maximum_iterations,  &
            & lrpar, rpar, lipar, ipar, ierr)
        if (ierr /= 0) then
            print *,'unable to bracket root: ierr = ', ierr
            Pphoto = fallback_Pphoto
            return
        end if
        
     	! set iteration count, tolerances
        maximum_iterations = default_maximum_iterations_photosphere
        eps_rho = default_tolerance_photosphere_rho
        eps_ph = default_tolerance_photosphere_condition

		root_ph = safe_root_with_initial_guess(photosphere,rho_guess,rho1,rho3,y1,y3, &
            &   maximum_iterations,eps_rho,eps_ph,lrpar,rpar,lipar,ipar,ierr)

        if (ierr /= 0) then
            write(*,*) 'unable to converge', rho_guess, rho1, rho3, y1, y3
            return
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
       real(dp) :: gravity, tau_ph, Teff, P, kappa, Gamma, eta, Xsum, chi
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
           ierr = negative_photosphere_density
           return
       end if
       
       gravity = rpar(igrav)
       tau_ph = rpar(itau)
       Teff = rpar(iTeff)
       eos_handle = ipar(ihandle)
	   chem_ids = ipar(ichem_id)
       chi = use_default_nuclear_size
       ! single species only
       Y(1) = 1.0/nuclib% A(chem_ids(1))

       call compute_composition_moments(1,chem_ids,Y,ionic,Xsum,ncharged, charged_ids, Yion,  &
   			& exclude_neutrons = .TRUE.)
       call eval_crust_eos(eos_handle,rho,Teff,ionic,ncharged,charged_ids,Yion, &
       		&   res,phase,chi)
       
       P = exp(res(i_lnP))
       rpar(iPph) = P
       eta = res(i_Theta) !1.0/TpT
       Gamma = res(i_Gamma)
       call get_thermal_conductivity(rho,Teff,chi, &
           & Gamma,eta,ionic,K,which_components=cond_use_only_kap)
       kappa = 4.0*onethird*arad*clight*Teff**3/rho/K% kap
       rpar(iKph) = kappa
	   dfdrho = 0.0
       
       photosphere = P - 2.0*onethird*gravity/kappa
    end function photosphere

end module bc09
