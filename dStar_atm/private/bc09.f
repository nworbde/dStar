module bc09
    use constants_def, only: dp
    
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
    integer, parameter :: ipres = iTeff + 1
    integer, parameter :: itemp = ipres + 1
    integer, parameter :: irho = itemp + 1
    integer, parameter :: ikappa = irho + 1
    integer, parameter :: number_base_rpar = ikappa
    ! further storage can be tacked on at end of rpar for Yion's
    
    integer, parameter :: ieoshandle = 1
    integer, parameter :: icondhandle = ieoshandle + 1
    integer, parameter :: iNcharged = icondhandle + 1
    integer, parameter :: number_base_ipar = iNcharged
    ! further storage can be tacked on at end of ipar for charged id's
    
    ! error codes
    integer, parameter :: negative_photosphere_gas_pressure = -2
    integer, parameter :: negative_photosphere_density = -3
    integer, parameter :: bad_composition = -4
        
    ! for debugging
    logical, parameter :: dbg = .FALSE.
    
contains
    
    subroutine do_get_bc09_Teff(grav, Plight, Pb, lgTb, lgTeff, lgflux, ierr)
        use constants_def
        use dStar_eos_lib
        use conductivity_lib
        
        real(dp), intent(in) :: grav    ! surface gravity, in the local frame
        real(dp), intent(in) :: Plight  ! pressure at which layer of light elements terminates
        real(dp), intent(in) :: Pb      ! pressure at base of layer
        real(dp), intent(out), dimension(:) :: lgTb ! temperature at a base column
        real(dp), intent(in), dimension(:) :: lgTeff    ! effective temperature
        real(dp), intent(out), dimension(:) :: lgflux   ! flux
        integer, intent(out) :: ierr
        integer :: i, n
        integer :: eos_handle, cond_handle
        real(dp) :: lgyb, lgy_light, rho_ph
        
        ! constants, nucchem, and dStar_eos must be initialized
        eos_handle = alloc_dStar_eos_handle(ierr)
        cond_handle = alloc_conductivity_handle(ierr)
        call conductivity_set_controls(cond_handle, &
        &   include_electrons=.TRUE., &
        &   include_photons=.TRUE., &
        &   include_neutrons=.FALSE., &
        &   include_superfluid_phonons=.FALSE.)
        
        lgyb = log10(Pb/grav)
        lgy_light = log10(Plight/grav)
        
        rho_ph = -1.0_dp
        ! start at the highest temperature and work down; array is assumed to be in ascending order
        n = size(lgTeff)
        do i = n,1,-1
            call do_integrate_bc09_atm( &
            &   grav,lgyb,lgy_light,lgTeff(i),lgTb(i),rho_ph, &
            &   eos_handle,cond_handle,ierr)
            if (ierr < 0) then
                print *,'error in atmosphere integration with lgTeff = ', &
                &    lgTeff(i)
                return
            end if
            lgflux(i) = 4.0_dp*lgTeff(i)+log10(sigma_SB)
        end do
        
        call free_conductivity_handle(cond_handle)
        call free_dStar_eos_handle(eos_handle)
    end subroutine do_get_bc09_Teff
    
    subroutine do_integrate_bc09_atm( &
    &   grav,lgyb,lgy_light,lgTeff,lgTb,rho,eos_handle,cond_handle,ierr)
        use iso_fortran_env, only: error_unit
        use constants_def
        use nucchem_def
        use nucchem_lib
        use num_lib
        
        real(dp), intent(in) :: grav    ! cm/s**2
        real(dp), intent(in) :: lgyb    ! log_10(g/cm**2); base of atmosphere
        real(dp), intent(in) :: lgy_light   ! log_10(g/cm**2); light/heavy transition
        real(dp), intent(in) :: lgTeff  ! log_10(K)
        real(dp), intent(out) :: lgTb   ! log_10(K)
        real(dp), intent(inout) :: rho  ! set < 0 to compute a guess, on output, contains rho_photosphere
        integer, intent(in) :: eos_handle, cond_handle
        integer, intent(out) :: ierr
        ! composition
        integer, parameter :: number_species = 2
        integer, dimension(number_species) :: charged_ids, chem_ids
        real(dp), dimension(number_species) :: Y, Yion
        real(dp) :: Xsum
        integer :: ncharged
        ! data arrays
        integer :: lrpar, lipar
        integer, pointer :: ipar(:) => null() ! (lipar)
        real(dp), pointer :: rpar(:) => null()  ! (lrpar)
        real(dp) :: Teff, tau, P, kappa
        type(composition_info_type) :: ionic
        ! integration
        real(dp) :: lnP, lnPend, h, max_step_size, atol(1), rtol(1)
        integer :: iout, itol, lwork, liwork, lout, idid, max_steps
        real(dp), pointer, dimension(:) :: lnT4=>null()
        integer, pointer, dimension(:) :: iwork=>null()
        real(dp), pointer, dimension(:) :: work=>null()
        
        ierr = 0
        ! composition is a He/Fe mix
        chem_ids = [ get_nuclide_index('he4'), get_nuclide_index('fe56') ]
        if (any(chem_ids == nuclide_not_found)) then    ! this is a fatal error
            ierr = nuclide_not_found
            return
        endif
        
        ! set size of data structure        
        lrpar = number_base_rpar + number_species
        lipar = number_base_ipar + number_species
        allocate(ipar(lipar), rpar(lrpar))
        
        Teff = 10.0**lgTeff
        tau = twothird
        
        ! set photosphere to be pure He and compute moments
        Y = [1.0_dp, 0.0_dp]/nuclib% A(chem_ids)
        call compute_composition_moments(number_species,chem_ids,Y,ionic,Xsum, &
        &   ncharged, charged_ids, Yion, exclude_neutrons = .TRUE.)

        ! sanity check on composition
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
        
        ipar(ieoshandle) = eos_handle
        ipar(icondhandle) = cond_handle
        ipar(iNcharged) = ncharged
        ipar(number_base_ipar+1:number_base_ipar+ncharged) = charged_ids(:)
        
        if (dbg) write (error_unit,*) 'finding photosphere...'
        call find_photospheric_pressure(Teff,grav,tau,rho,P,kappa, &
            &   lrpar, rpar, lipar, ipar, ierr)
        if (ierr /= 0) return

        if (dbg) then
            write (error_unit,*) 'done'
            write (error_unit,'(2(a,es11.4))') 'rho = ',rho,'; kappa = ',kappa
        end if
        ! anchor the eos information
        rpar(irho) = rho
        rpar(itemp) = Teff
        rpar(ipres) = P
        rpar(ikappa) = kappa

        ! start the integration
        if (dbg) write (error_unit,*) 'starting integration'
        call dop853_work_sizes(1,0,liwork,lwork)
        allocate(lnT4(1),iwork(liwork),work(lwork))
        
        ! integration variables
        lnP = log(P)
        lnT4(1) = 0.0   ! 4.0*ln(T/Teff)
        
        ! first layer; we shall integrate over a minimum of 1 decade in pressure
        lnPend = max(ln10*(lgy_light + log10(grav)), lnP+ln10)
        ! but the light element layer cannot exceed the total thickness
        lnPend = min(lnPend, ln10*(lgyb + log10(grav)))
        
        h = 0.001_dp
        max_step_size = 0.0_dp
        max_steps = 10000
        rtol(:) = 1.0e-6_dp
        atol(:) = 1.0e-6_dp
        itol = 0
        iout = 0
        lout = error_unit
        iwork(:) = 0
        work(:) = 0.0
        call dop853(1,deriv,lnP,lnT4,lnPend,h,max_step_size,max_steps, &
        &   rtol, atol, itol, null_solout, iout, work, lwork, iwork, liwork, &
        &   lrpar, rpar, lipar, ipar, lout, idid)
        
        if (idid < 0) then
            write (error_unit,*) 'error in integration'
            write (error_unit,'(2(a,f8.4))') 'lnP = ',lnP,'; lnT4 = ',lnT4(1)
            ierr = idid 
            return
        end if

        if (dbg) write (error_unit,*) 'integrating over heavy layer'
        ! now integrate over the heavy layer
        Y = [0.0_dp, 1.0_dp]/nuclib% A(chem_ids)
        call compute_composition_moments(number_species,chem_ids,Y,ionic,Xsum, &
        &   ncharged, charged_ids, Yion, exclude_neutrons = .TRUE.)

        ! sanity check on composition
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
        rpar(number_base_rpar + 1:number_base_rpar+ncharged) = Yion(:)        
        ipar(number_base_ipar+1:number_base_ipar+ncharged) = charged_ids(:)
 
        ! reset end point
        lnPend = ln10*(lgyb + log10(grav))
        
        if (lnP < lnPend) then  ! integrate over the layer
            ! reset the step size since we have a new composition
            h = 0.001_dp
            iwork(:) = 0
            work(:) = 0.0_dp
            call dop853(1,deriv,lnP,lnT4,lnPend,h,max_step_size,max_steps, &
            &   rtol, atol, itol, null_solout, iout, work, lwork, iwork, liwork, &
            &   lrpar, rpar, lipar, ipar, lout, idid)
        
            if (idid < 0) then
                write (error_unit,*) 'error in integration'
                ierr = idid
                return
            end if
        end if
        lgTb = 0.25_dp*lnT4(1)/ln10 + lgTeff
                
        deallocate(ipar, rpar, iwork, work, lnT4)
    end subroutine do_integrate_bc09_atm
      
    subroutine find_photospheric_pressure(Teff,grav,tau,rho_ph,P_ph,kappa, &
        &   lrpar,rpar,lipar,ipar,ierr)
        use iso_fortran_env, only: error_unit
        use constants_def
        use nucchem_def
        use nucchem_lib
        use num_lib
        integer, parameter :: default_maximum_iterations_photosphere = 20
        real(dp), parameter :: default_tolerance_photosphere_lnrho = 1.0e-6_dp
        real(dp), parameter :: default_tolerance_photosphere_condition = 1.0e-8_dp
        real(dp), intent(in) :: Teff    ! K
        real(dp), intent(in) :: grav    ! cm/s**2
        real(dp), intent(in) :: tau     ! may need to adjust to something other than 2.0/3.0, 
            ! especially at higher temperatures
        real(dp), intent(inout) :: rho_ph   ! on input set <= 0 to have routine generate guess; on 
        !   output, it contains the value of the photospheric density
        real(dp), intent(out) :: P_ph   ! photospheric pressure, cgs units
        real(dp), intent(out) :: kappa  ! opacity at photosphere
        integer, intent(in) :: lrpar, lipar
        integer, intent(inout), pointer :: ipar(:) ! (lipar)
        real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
        integer, intent(out) :: ierr
        real(dp) :: lnrho_ph,lnrho
        integer :: i, maximum_iterations
        real(dp) :: sigma_Th
        real(dp) :: lnrho_guess, kappa_Th, Pgas, fallback_Pphoto
        real(dp) :: lnrho1, lnrho3, dlnrho, ph1, ph3, eps_lnrho, eps_ph
        
        ierr = 0

        ! set iteration count, tolerances
        maximum_iterations = default_maximum_iterations_photosphere
        eps_lnrho = default_tolerance_photosphere_lnrho
        eps_ph = default_tolerance_photosphere_condition
        fallback_Pphoto = 2.0_dp*onethird*arad*Teff**4

        ! scale tolerance to a thomson scattering atmosphere
        sigma_Th = 8.0_dp*onethird*pi*(electroncharge**2/Melectron/clight2)**2
        kappa_Th = sigma_Th*avogadro*rpar(icomp_Ye)
        eps_ph = default_tolerance_photosphere_condition*tau*grav/kappa_Th
        
        ! use initial guess with ideal gas pressure and thomson scattering
        if (rho_ph < 0.0_dp) then
            Pgas = tau*grav/kappa_Th - onethird*arad*Teff**4
            if (Pgas < 0.0_dp) then
                write(error_unit,*) 'negative photosphere gas pressure'
                ierr = negative_photosphere_gas_pressure
                P_ph = fallback_Pphoto
                kappa = kappa_Th
                return
            end if
            lnrho_guess = log(Pgas*amu*rpar(icomp_A)/(rpar(icomp_Z)+1.0)/boltzmann/Teff)
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
            kappa = rpar(ikappa)
            return
        end if

        lnrho_ph = safe_root_with_initial_guess(photosphere,lnrho_guess,lnrho1,lnrho3,ph1,ph3, &
            &   maximum_iterations,eps_lnrho,eps_ph,lrpar,rpar,lipar,ipar,ierr)

        if (ierr /= 0) then
            write(error_unit,*) 'unable to converge on photospheric density: ',lnrho_guess,lnrho_ph
            P_ph = fallback_Pphoto
            rho_ph = exp(lnrho_guess)
            kappa = rpar(ikappa)
            return
        end if
        P_ph = rpar(ipres)
        rho_ph = exp(lnrho_ph)
        kappa = rpar(ikappa)
    end subroutine find_photospheric_pressure    

    real(dp) function photosphere(lnrho, dfdlnrho, lrpar, rpar, lipar, ipar, ierr)
       ! returns with ierr = 0 if was able to evaluate f and df/dx at x
       ! if df/dx not available, it is okay to set it to 0
       use constants_def
       use superfluid_def, only: max_number_sf_types, neutron_1S0
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
       real(dp) :: rho, gravity, tau_ph, Teff, P, kappa, Gamma, eta, mu_e, Xsum, chi
       integer :: eos_handle, cond_handle, ncharged, phase
       type(composition_info_type) :: ionic
       type(conductivity_components) :: K
       real(dp), dimension(num_dStar_eos_results) :: res
       integer, pointer, dimension(:) :: charged_ids=>null()
       real(dp), pointer, dimension(:) :: Yion=>null()
       real(dp), dimension(max_number_sf_types) :: Tcs
           
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
       
       eos_handle = ipar(ieoshandle)
       cond_handle = ipar(icondhandle)
       ncharged = ipar(iNcharged)
       charged_ids(1:ncharged)=>ipar(number_base_ipar+1:number_base_ipar+ncharged)
       Yion(1:ncharged)=>rpar(number_base_rpar+1:number_base_rpar+ncharged)
       chi = use_default_nuclear_size

       Tcs = 0.0_dp
       call eval_crust_eos(eos_handle,rho,Teff,ionic,ncharged,charged_ids, &
        &   Yion, Tcs, res,phase,chi)
       
       P = exp(res(i_lnP))
       rpar(ipres) = P
       eta = res(i_Theta) !1.0/TpT
       mu_e = res(i_mu_e)
       Gamma = res(i_Gamma)
       call get_thermal_conductivity(cond_handle,rho,Teff,chi, &
           & Gamma,eta,mu_e,ionic,Tcs(neutron_1S0), &
           & K)
       kappa = 4.0*onethird*arad*clight*Teff**3/rho/K% total
       rpar(ikappa) = kappa
       rpar(ipres) = P
       dfdlnrho = 0.0
       photosphere = P - tau_ph*gravity/kappa
    end function photosphere
    
    subroutine deriv(n, lnP, h, lnT4, dlnT4dlnP, lrpar, rpar, lipar, ipar, ierr)
       use constants_def
       use nucchem_def
       use nucchem_lib
       use dStar_eos_lib
       use conductivity_lib
       use utils_lib

       integer, intent(in) :: n, lrpar, lipar
       real(dp), intent(in) :: lnP, h
       real(dp), intent(inout) :: lnT4(:)   ! 4.0*ln(T/Teff)
       real(dp), intent(inout) :: dlnT4dlnP(:)
       integer, intent(inout), pointer :: ipar(:) ! (lipar)
       real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
       integer, intent(out) :: ierr ! nonzero means retry with smaller timestep.
       real(dp) :: Teff, grav, rho, P, T, kappa, del_ad
       
       ierr = 0
       ! unpack the arguments
       Teff = rpar(iTeff)
       grav = rpar(igrav)
       rho = rpar(irho)

       P = exp(lnP)
       T = Teff*exp(0.25*lnT4(1))
       
       call get_coefficients(P,T,rho,lrpar,rpar,lipar,ipar,kappa,del_ad,ierr)
       if (ierr /= 0) then
           print *,'unable to get coefficients'
           return
       end if
       
       dlnT4dlnP(1) = 0.75_dp*kappa*P/grav/exp(lnT4(1))
       dlnT4dlnP(1)  = min(dlnT4dlnP(1), 4.0*del_ad)       

       ! store the new anchor point
       rpar(irho) = rho

    end subroutine deriv

    subroutine get_coefficients(P,T,rho,lrpar,rpar,lipar,ipar,kappa,del_ad,ierr)
        use constants_def
        use superfluid_def, only: max_number_sf_types, neutron_1S0
        use nucchem_def, only: composition_info_type
        use dStar_eos_lib
        use conductivity_lib
        
        real(dp), intent(in) :: P,T
        real(dp), intent(inout) :: rho  ! on input, holds the guess for rho
        integer, intent(in) :: lrpar,lipar
        real(dp), dimension(:), pointer :: rpar
        integer, dimension(:), pointer :: ipar
        real(dp), intent(out) :: kappa
        real(dp), intent(out) :: del_ad
        integer, intent(out) :: ierr
        type(composition_info_type) :: ionic
        integer :: ncharged, eos_handle, cond_handle
        integer, dimension(:), pointer :: charged_ids=>null()
        real(dp), dimension(:), pointer :: Yion=>null()
        real(dp) :: lnrho,lnrho_guess, Gamma, eta, mu_e, chi
        real(dp), dimension(max_number_sf_types) :: Tcs
        real(dp), dimension(num_dStar_eos_results) :: res
        integer :: phase
        type(conductivity_components) :: K
        
        lnrho_guess = log(rho)
        lnrho = get_lnrho_from_PT(P,T,lnrho_guess,lrpar,rpar,lipar,ipar,ierr)
        if (ierr /= 0) then
            print *,'unable to get density'
            return
        end if
        
        rho = exp(lnrho)
        chi = use_default_nuclear_size
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
        eos_handle = ipar(ieoshandle)
        cond_handle = ipar(icondhandle)
        ncharged = ipar(iNcharged)
        charged_ids(1:ncharged) => ipar(number_base_ipar+1:number_base_ipar+ncharged)
        Yion(1:ncharged) => rpar(number_base_rpar+1:number_base_rpar+ncharged)
        Tcs = 0.0_dp
        call eval_crust_eos(eos_handle,rho,T,ionic,ncharged,charged_ids,Yion, &
                &   Tcs, res,phase,chi)
                
        ! check for good eos
        if (abs(exp(res(i_lnP))-P) > 1.0e-3_dp*P) then
            print *, 'bad rootfind: ',rho,T,P,exp(res(i_lnP)),ionic% Z, ionic %A
            stop
        end if
        
        eta = res(i_Theta) !1.0/TpT
        Gamma = res(i_Gamma)
        mu_e = res(i_mu_e)
        del_ad = res(i_grad_ad)
        call get_thermal_conductivity(cond_handle,rho,T,chi, &
            & Gamma,eta,mu_e,ionic,Tcs(neutron_1S0), &
            & K)
        kappa = 4.0*onethird*arad*clight*T**3/rho/K% total
        
!         if (dbg) then
!             print *,rho,eta,Gamma,mu_e,del_ad,kappa
!         end if
    end subroutine get_coefficients

    function get_lnrho_from_PT(P,T,lnrho_guess,lrpar,rpar,lipar,ipar,ierr) result(lnrho)
        use iso_fortran_env, only: error_unit
        use constants_def
        use nucchem_def, only: composition_info_type
        use dStar_eos_lib
        use num_lib, only: look_for_brackets, safe_root_with_initial_guess
        real(dp), intent(in) :: P
        real(dp), intent(in) :: T
        real(dp), intent(in) :: lnrho_guess
        integer, intent(in) :: lrpar,lipar
        real(dp), dimension(:), pointer :: rpar
        integer, dimension(:), pointer :: ipar
        integer, intent(out) :: ierr
        real(dp) :: lnrho
        real(dp) :: dlnrho  ! increment for searching for brackets
        real(dp) :: lnrho1, lnrho3, p1, p3, dlnPdlnrho
        integer, parameter :: default_maximum_iterations_density = 20
        real(dp), parameter :: default_tolerance = 1.0e-12_dp
        integer :: maximum_iterations
        real(dp) :: eps_rho,eps_p
        
        ierr = 0
                    
        maximum_iterations = default_maximum_iterations_density
        eps_rho = default_tolerance
        eps_p = default_tolerance
        
        ! store temp, pressure for root find
        rpar(itemp) = T
        rpar(ipres) = P
        
        ! get brackets for root find
        dlnrho = 0.5_dp
        call look_for_brackets(lnrho_guess, dlnrho, lnrho1, lnrho3,  &
        &   eval_pressure, p1, p3, maximum_iterations,  &
        &   lrpar, rpar, lipar, ipar, ierr)
        if (ierr /= 0) then
            write(error_unit,*) 'unable to bracket density: ierr = ', ierr
            return
        end if
        
        lnrho = safe_root_with_initial_guess(eval_pressure,lnrho_guess,lnrho1,lnrho3,p1,p3, &
            &   maximum_iterations,eps_rho,eps_p,lrpar,rpar,lipar,ipar,ierr)

        if (ierr /= 0) then
            write(error_unit,'(a,2(es13.4))') 'unable to converge to density: ',lnrho_guess,lnrho
            return
        end if

!         print *,'eval p =', eval_pressure(lnrho,dlnPdlnrho,lrpar,rpar,lipar,ipar,ierr)
        if (ierr /=0 ) then
            print *,'problem in eval_pressure'
            return
        end if
        
    end function get_lnrho_from_PT

    real(dp) function eval_pressure(lnrho, dlnPdlnrho, lrpar, rpar, lipar, ipar, ierr)
       ! returns with ierr = 0 if was able to evaluate lnP and dlnP/dlnrho at rho
       use constants_def
       use superfluid_def, only: max_number_sf_types
       use nucchem_def
       use nucchem_lib
       use dStar_eos_lib
       use utils_lib

       integer, intent(in) :: lrpar, lipar
       real(dp), intent(in) :: lnrho
       real(dp), intent(out) :: dlnPdlnrho
       integer, intent(inout), pointer :: ipar(:) ! (lipar)
       real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
       integer, intent(out) :: ierr
       real(dp) :: rho, T, chi, Pwant
       integer :: eos_handle, ncharged, phase
       type(composition_info_type) :: ionic
       real(dp), dimension(max_number_sf_types) :: Tcs
       real(dp), dimension(num_dStar_eos_results) :: res
       integer, pointer, dimension(:) :: charged_ids=>null()
       real(dp), pointer, dimension(:) :: Yion=>null()

       ierr = 0
       rho = exp(lnrho)
       T = rpar(itemp)
       Pwant = rpar(ipres)
       eos_handle = ipar(ieoshandle)
       ncharged = ipar(iNcharged)
       charged_ids(1:ncharged)=>ipar(number_base_ipar+1:number_base_ipar+ncharged)
       Yion(1:ncharged)=>rpar(number_base_rpar+1:number_base_rpar+ncharged)
       chi = use_default_nuclear_size       
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

       if (sum(Yion*nuclib% A(charged_ids)) - 1.0_dp > 4.0*epsilon(1.0_dp)) then
           print *,'bad composition:', Yion(1:ncharged), charged_ids(1:ncharged)
           stop
       end if
       
       Tcs = 0.0_dp
       call eval_crust_eos(eos_handle,rho,T,ionic,ncharged,charged_ids,Yion, &
               &   Tcs, res,phase,chi)
       if (is_bad_num(res(i_lnP))) then
           ierr = -9
           return
       end if
       dlnPdlnrho = res(i_chiRho)
       eval_pressure = res(i_lnP) - log(Pwant)
    end function eval_pressure

end module bc09
