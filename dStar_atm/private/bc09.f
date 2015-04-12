module bc09
    
    ! data structure for the parameter arrays
    ! moments of the composition
    integer, parameter :: icomp_A = 1
    integer, parameter :: icomp_Z = icomp_A + 1
    integer, parameter :: icomp_Z53 = icomp_Z + 1
    integer, parameter :: icomp_Z2 = icomp_Z53 + 1
    integer, parameter :: icomp_Z73 = icomp_Z2 + 1
    integer, parameter :: icomp_Z52 = icomp_Z73 + 1
    integer, parameter :: icomp_ZZ1_32 = icomp_Z52 + 1
    integer, parameter :: icomp_Z2XoA2 = icomp_ZZ1_32 + 1
    integer, parameter :: icomp_Ye = icomp_Z2XoA2 + 1
    integer, parameter :: icomp_Yn = icomp_Ye + 1
    integer, parameter :: icomp_Q = icomp_Yn + 1
    integer, parameter :: number_comp_rpar = icomp_Q
        
    integer, parameter :: igrav = number_comp_rpar + 1
    integer, parameter :: itau = igrav + 1
    integer, parameter :: iTeff = itau + 1
    integer, parameter :: iPph = iTeff + 1
    integer, parameter :: iKph = iPph + 1
    integer, parameter :: iChi_rho = iKph + 1
    integer, parameter :: iChi_T = iChi_rho + 1
    integer, parameter :: del_ad = iChi_T + 1
    integer, parameter :: number_base_rpar = del_ad
    ! further storage can be tacked on at end of rpar for Yion's
    
    integer, parameter :: ihandle = 1
    integer, parameter :: iNcharged = ihandle + 1
    integer, parameter :: number_base_ipar = iNcharged
    ! further storage can be tacked on at end of ipar for charged id's
    
    ! error codes

    integer, parameter :: negative_photosphere_gas_pressure = -2
    integer, parameter :: negative_photosphere_density = -3
    integer, parameter :: bad_composition = -4
contains
    
	subroutine do_get_bc09_Teff(grav, Plight, Tb, Teff, flux)
		use constants_def
		real(dp), intent(in) :: grav	! surface gravity, in the local frame
		real(dp), intent(in) :: Plight	! pressure at which layer of light elements terminates
		real(dp), intent(in), dimension(:) :: Tb	! temperature at a base column
		real(dp), intent(out), dimension(:) :: Teff, flux	! effective temperature and flux
		real(dp) :: eta, g14
        integer ::  size_tab ! = 4*size(Tb)
        real(dp), dimension(:), allocatable :: tabTb9, tabTeff, tabTeff6_4, tabRhoph, tabPph
        integer :: i
        
        ! make a very dense table of Tb(Teff); then interpolate to get Teff(Tb)        
        size_tab = 4*size(Tb)
        allocate(tabTb9(size_tab),tabTeff(size_tab),tabTeff6_4(size_tab) &
        	&	tabRhoph(size_tab), tabPph(size_tab))

    	!tau = ???
		rho_ph = -1.0
        ! compute dense table; get Pph(Teff)
        do i = size_tab,1,-1
       	 	lnTeff = log10(5.0e5) + (i-1)/real(size_tab-1)
        	Teff = 10.0_dp**lnTeff
        	tabTeff(i) = Teff
        	tabTeff6_4(i) = (Teff/1.0d6)**4.0
        	call find_photospheric_pressure(Teff,grav,tau,rho_ph,P_ph,kappa,eos_handle,ierr)
        	if (ierr /= 0) then
           		print *,'error: ierr = ',ierr
            	rho_ph = -1.0
            	cycle
       		end if
		    tabRhoph(i) = rho_ph
		    tabPph(i) = P_ph
    	end do        
               
        ! interpolate from dense table to get finished product
        
        deallocate(tabTb9,tabTeff,tabTeff6_4)
    end subroutine do_get_bc09_Teff
    
    subroutine find_photospheric_pressure(Teff,grav,tau,rho_ph,P_ph,kappa,eos_handle,ierr)
        use iso_fortran_env, only: error_unit
        use constants_def
        use nucchem_def
    	use nucchem_lib
    	use num_lib
        integer, parameter :: default_maximum_iterations_photosphere = 20
        real(dp), parameter :: default_tolerance_photosphere_lnrho = 1.0e-6_dp
        real(dp), parameter :: default_tolerance_photosphere_condition = 1.0e-8_dp
        integer, parameter :: number_species = 2
        real(dp), intent(in) :: Teff,grav,tau
        real(dp), intent(inout) :: rho_ph   ! on input set <= 0 to have routine generate guess; on 
        !   output, it contains the value of the photospheric density
        real(dp), intent(out) :: P_ph,kappa
        integer, intent(in) :: eos_handle
        integer, intent(out) :: ierr
        real(dp) :: lnrho_ph,lnrho
        integer, pointer :: ipar(:) => null() ! (lipar)
        real(dp), pointer :: rpar(:) => null()  ! (lrpar)
        integer, dimension(number_species) :: charged_ids, chem_ids
        real(dp), dimension(number_species) :: Y, Yion
        integer :: i, lrpar, lipar, maximum_iterations, ncharged
        real(dp) :: sigma_Th, Xsum
        real(dp) :: lnrho_guess, kappa_Th, Pgas, fallback_Pphoto
        real(dp) :: lnrho1, lnrho3, dlnrho, ph1, ph3, eps_lnrho, eps_ph
        type(composition_info_type) :: ionic
        
        ierr = 0
        
        ! composition is a He/Fe mix
        ! set size of data structure for rpar, ipar
        
        lrpar = number_base_rpar + number_species
        lipar = number_base_ipar + number_species
        allocate(ipar(lipar), rpar(lrpar))

        chem_ids = [ get_nuclide_index('he4'), get_nuclide_index('fe56') ]
        if (any(chem_ids == nuclide_not_found)) then    ! this is a fatal error
            ierr = nuclide_not_found
            return
        endif
        
        ! set photosphere to be pure He
        Y = [1.0_dp, 0.0_dp]
        Y = Y/nuclib% A(chem_ids)

        call compute_composition_moments(number_species,chem_ids,Y,ionic,Xsum, &
        &   ncharged, charged_ids, Yion, exclude_neutrons = .TRUE.)

        if (Xsum - 1.0_dp > 2.0*epsilon(1.0_dp)) then
            ierr = bad_composition
            return
        end if
        
        ! stuff the parameter vectors
        rpar(icomp_A) = ionic% A
        rpar(icomp_Z) = ionic% Z
        rpar(icomp_Z53) = ionic% Z53
        rpar(icomp_Z2) = ionic% Z2
        rpar(icomp_Z73) = ionic% Z73
        rpar(icomp_Z52) = ionic% Z52
        rpar(icomp_ZZ1_32) = ionic% ZZ1_32
        rpar(icomp_Z2XoA2) = ionic% Z2XoA2
        rpar(icomp_Ye) = ionic% Ye
        rpar(icomp_Yn) = ionic% Yn
        rpar(icomp_Q) = ionic% Q
        rpar(igrav) = grav
        rpar(itau) = tau
        rpar(iTeff) = Teff
        rpar(number_base_rpar + 1:number_base_rpar+ncharged) = Yion(:)
        
        ipar(ihandle) = eos_handle
        ipar(iNcharged) = ncharged
        ipar(number_base_ipar+1:number_base_ipar+ncharged) = charged_ids

     	! set iteration count, tolerances
        maximum_iterations = default_maximum_iterations_photosphere
        eps_lnrho = default_tolerance_photosphere_lnrho
        eps_ph = default_tolerance_photosphere_condition
        fallback_Pphoto = 2.0_dp*onethird*arad*Teff**4   

        ! scale tolerance to a thomson scaterring atmosphere
        sigma_Th = 8.0_dp*onethird*pi*(electroncharge**2/Melectron/clight2)**2
        kappa_Th = sigma_Th*avogadro*ionic% Ye
        eps_ph = default_tolerance_photosphere_condition*tau*grav/kappa_Th
        
        ! use initial guess with ideal gas pressure and thomson scattering
        if (rho_ph < 0.0_dp) then
            Pgas = tau*grav/kappa_Th - onethird*arad*Teff**4
            if (Pgas < 0.0_dp) then
                ierr = negative_photosphere_gas_pressure
                P_ph = fallback_Pphoto
                kappa = kappa_Th
                return
            end if
            lnrho_guess = log(Pgas*amu*ionic% A/(ionic% Z+1.0)/boltzmann/Teff)
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
            kappa = rpar(iKph)
            return
        end if
        
		lnrho_ph = safe_root_with_initial_guess(photosphere,lnrho_guess,lnrho1,lnrho3,ph1,ph3, &
            &   maximum_iterations,eps_lnrho,eps_ph,lrpar,rpar,lipar,ipar,ierr)

        if (ierr /= 0) then
            write(error_unit,*) 'unable to converge on photospheric density: ',lnrho_guess,lnrho_ph
            P_ph = fallback_Pphoto
            rho_ph = exp(lnrho_guess)
            kappa = rpar(iKph)
            return
        end if
        P_ph = rpar(iPph)
        rho_ph = exp(lnrho_ph)
        kappa = rpar(iKph)
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
       integer :: eos_handle, ncharged, phase
       type(composition_info_type) :: ionic
       type(conductivity_components) :: K
       real(dp), dimension(num_dStar_eos_results) :: res
	   integer, pointer, dimension(:) :: charged_ids=>null()
       real(dp), pointer, dimension(:) :: Yion=>null()
       
       ierr = 0
       
       ! unpack the arguments       
       rho = exp(lnrho)
       
       ionic = composition_info_type(A = rpar(icomp_A), &
       &    Z = rpar(icomp_Z), &
       &    Z53 = rpar(icomp_Z53), &
       &    Z2 = rpar(icomp_Z2), &
       &    Z73 = rpar(icomp_Z73), &
       &    Z52 = rpar(icomp_Z52), &
       &    ZZ1_32 = rpar(icomp_ZZ1_32), &
       &    Z2XoA2 = rpar(icomp_Z2XoA2), &
       &    Ye = rpar(icomp_Ye), &
       &    Yn = rpar(icomp_Yn), &
       &    Q = rpar(icomp_Q) )
       
       gravity = rpar(igrav)
       tau_ph = rpar(itau)
       Teff = rpar(iTeff)
       
       eos_handle = ipar(ihandle)
       ncharged = ipar(iNcharged)
       charged_ids=>ipar(number_base_ipar+1:number_base_ipar+ncharged)
       Yion=>rpar(number_base_rpar+1:number_base_rpar+ncharged)
       chi = use_default_nuclear_size

       call eval_crust_eos(eos_handle,rho,Teff,ionic,ncharged,charged_ids,Yion, &
       		&   res,phase,chi)
       
       P = exp(res(i_lnP))
       rpar(iPph) = P
       eta = res(i_Theta) !1.0/TpT
       Gamma = res(i_Gamma)
       call get_thermal_conductivity(rho,Teff,chi, &
           & Gamma,eta,ionic,K,which_components=cond_exclude_sf) !cond_use_only_kap)
       kappa = 4.0*onethird*arad*clight*Teff**3/rho/K% kap
       rpar(iKph) = kappa
	   dfdlnrho = 0.0
       
       photosphere = P - 2.0*onethird*gravity/kappa
    end function photosphere

end module bc09
