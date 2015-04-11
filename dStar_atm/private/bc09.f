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
    
    subroutine find_photospheric_pressure(Teff,grav,tau,rho_ph,P_ph,eos_handle,ierr)
        use iso_fortran_env, only: error_unit
        use constants_def
        use nucchem_def, only : nuclide_not_found, nuclib
    	use nucchem_lib, only : get_nuclide_index
    	use num_lib
        integer, parameter :: default_maximum_iterations_photosphere = 20
        real(dp), parameter :: default_tolerance_photosphere_lnrho = 1.0e-6_dp
        real(dp), parameter :: default_tolerance_photosphere_condition = 1.0e-6_dp
        real(dp), intent(in) :: Teff,grav,tau
        real(dp), intent(inout) :: rho_ph   ! on input set <= 0 to have routine generate guess; on 
        !   output, it contains the value of the photospheric density
        real(dp), intent(out) :: P_ph
        integer, intent(in) :: eos_handle
        integer, intent(out) :: ierr
        real(dp) :: lnrho_ph,lnrho
        integer, pointer :: ipar(:) => null() ! (lipar)
        real(dp), pointer :: rpar(:) => null()  ! (lrpar)
        integer :: i, lrpar, lipar
        real(dp) :: sigma_Th
        real(dp) :: fallback_Pphoto
        real(dp) :: lnrho_guess, kappa_Th, A, Z, Pgas
        real(dp) :: lnrho1, lnrho3, dlnrho, ph1, ph3
        real(dp) :: eps_lnrho, eps_ph
        integer :: maximum_iterations
        
        ierr = 0
        allocate(ipar(number_photosphere_ipar), rpar(number_photosphere_rpar))

        lrpar = number_photosphere_rpar
        lipar = number_photosphere_ipar
        rpar(igrav) = grav
        rpar(itau) = tau
        rpar(iTeff) = Teff
        ipar(ihandle) = eos_handle

        ! hardwire photosphere compostion to pure He
        ipar(ichem_id) = get_nuclide_index('he4')
        if (ipar(ichem_id) == nuclide_not_found) then
            ierr = nuclide_not_found
            return
        end if

     	! set iteration count, tolerances
        maximum_iterations = default_maximum_iterations_photosphere
        eps_lnrho = default_tolerance_photosphere_lnrho
        eps_ph = default_tolerance_photosphere_condition
        fallback_Pphoto = 2.0_dp*onethird*arad*Teff**4   
        
        ! use initial guess with ideal gas pressure and thomson scattering
        if (rho_ph < 0.0_dp) then
            A = nuclib% A(ipar(ichem_id))
            Z = nuclib% Z(ipar(ichem_id))
            sigma_Th = 8.0_dp*onethird*pi*(electroncharge**2/Melectron/clight2)**2
            kappa_Th = sigma_Th/amu* Z/A
            Pgas = tau*grav/kappa_th - onethird*arad*Teff**4
            if (Pgas < 0.0_dp) then
                ierr = negative_photosphere_gas_pressure
                P_ph = fallback_Pphoto
                return
            end if
            lnrho_guess = log(Pgas*amu*A/(Z+1.0)/(boltzmann*Teff))
        else
            lnrho_guess = log(rho_ph)
        end if

        ! get brackets for root find
        dlnrho = 0.1_dp
        call look_for_brackets(lnrho_guess, dlnrho, lnrho1, lnrho3, photosphere, ph1, ph3, &
             & maximum_iterations, lrpar, rpar, lipar, ipar, ierr)
        if (ierr /= 0) then
            write(error_unit,*) 'unable to bracket root: ierr = ', ierr
            P_ph = fallback_Pphoto
            rho_ph = exp(lnrho_guess)
            return
        end if
        
		lnrho_ph = safe_root_with_initial_guess(photosphere,lnrho_guess,lnrho1,lnrho3,ph1,ph3, &
            &   maximum_iterations,eps_lnrho,eps_ph,lrpar,rpar,lipar,ipar,ierr)

        if (ierr /= 0) then
            write(error_unit,*) 'unable to converge on photospheric density: ',lnrho_guess,lnrho_ph
            P_ph = fallback_Pphoto
            rho_ph = exp(lnrho_guess)
            return
        end if
        P_ph = rpar(iPph)
        rho_ph = exp(lnrho_ph)
        deallocate(ipar, rpar)
    end subroutine find_photospheric_pressure    

    real(dp) function photosphere(lnrho, dfdlnrho, lrpar, rpar, lipar, ipar, ierr)
       ! returns with ierr = 0 if was able to evaluate f and df/dx at x
       ! if df/dx not available, it is okay to set it to 0
       use constants_def
       use nucchem_def
       use nucchem_lib
       use dStar_eos_lib
       use conductivity_lib
       
       integer, intent(in) :: lrpar, lipar
       real(dp), intent(in) :: lnrho
       real(dp), intent(out) :: dfdlnrho
       integer, intent(inout), pointer :: ipar(:) ! (lipar)
       real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
       integer, intent(out) :: ierr
       real(dp) :: rho, gravity, tau_ph, Teff, P, kappa, Gamma, eta, Xsum, chi
       integer :: eos_handle
       type(composition_info_type) :: ionic
       type(conductivity_components) :: K
       real(dp), dimension(num_dStar_eos_results) :: res
       integer :: phase, ncharged
       integer, dimension(1) :: charged_ids, chem_ids
	   real(dp),dimension(1) :: Y,Yion
	   
       ierr = 0
       
       rho = exp(lnrho)
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
	   dfdlnrho = 0.0
       
       photosphere = P - 2.0*onethird*gravity/kappa
    end function photosphere

end module bc09
