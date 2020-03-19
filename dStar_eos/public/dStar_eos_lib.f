module dStar_eos_lib
    use dStar_eos_def
    
    contains
    subroutine dStar_eos_startup(datadir)
        use exceptions_lib
        use dStar_eos_private_def
        use helm_alloc
        character(len=*), intent(in) :: datadir
        integer, parameter :: imax = 261, jmax = 101  ! dimensions of our version of helm table
        integer :: ierr
        character(len=128) :: eos_datadir
        type(failure) :: startup_error=failure(scope='dStar_eos_startup')
        
        call dStar_eos_def_init
        call alloc_helm_table(eos_ht, imax, jmax, ierr)
        if (startup_error% raised(ierr,'unable to alloc helm table')) return
        
        eos_datadir = trim(datadir)//'/eos'
        call read_helm_table(eos_ht,eos_datadir,ierr)
        if (startup_error% raised(ierr,'unable to read helm table')) return
    end subroutine dStar_eos_startup
    
    subroutine dStar_eos_shutdown()
        use helm_alloc, only: free_helm_table
        call free_helm_table(eos_ht)
    end subroutine dStar_eos_shutdown

    function alloc_dStar_eos_handle(ierr)
        use, intrinsic :: iso_fortran_env, only: error_unit
        use exceptions_lib
        use dStar_eos_private_def
        integer, intent(out) :: ierr
        integer :: alloc_dStar_eos_handle
        type(failure) :: alloc_error=failure(scope='alloc_dStar_eos_handle')
        alloc_dStar_eos_handle = do_alloc_dStar_eos(ierr)
        if (alloc_error% raised(ierr, &
        &   trim(dstar_eos_private_def_errors(ierr)))) return
    end function alloc_dStar_eos_handle

    subroutine free_dStar_eos_handle(handle)
        use dStar_eos_private_def
        integer, intent(in) :: handle
        call do_free_dStar_eos(handle)
    end subroutine free_dStar_eos_handle

    subroutine dStar_eos_ptr(handle, rq, ierr)
        use, intrinsic :: iso_fortran_env, only: error_unit
        use exceptions_lib
        use dStar_eos_private_def
        integer, intent(in) :: handle
        type(dStar_eos_general_info), pointer :: rq
        integer, intent(out) :: ierr
        type(failure) :: ptr_error=failure(scope='dStar_eos_ptr')
        call get_dStar_eos_ptr(handle,rq,ierr)
        if (ptr_error% raised(ierr,trim(dstar_eos_private_def_errors(ierr)))) &
        &   return
    end subroutine dStar_eos_ptr
    
    subroutine dStar_eos_set_controls(handle,gamma_melt_pt,rsi_melt_pt,nuclide_abundance_threshold, &
        &       pasta_transition_in_fm3,cluster_transition_in_fm3,suppress_warnings)
        use, intrinsic :: iso_fortran_env, only: error_unit
        use exceptions_lib
        use dStar_eos_private_def, only : dStar_eos_general_info
        integer, intent(in) :: handle
        real(dp), intent(in), optional :: gamma_melt_pt, rsi_melt_pt, nuclide_abundance_threshold
        real(dp), intent(in), optional :: pasta_transition_in_fm3, cluster_transition_in_fm3
        logical, intent(in), optional :: suppress_warnings
        type(dStar_eos_general_info), pointer :: rq
        integer :: ierr
        type(assertion) :: handle_allocated= &
        &   assertion(scope='dStar_eos_set_controls', &
        &   message='unallocated handle passed to dStar_eos_set_controls')
        
        call dStar_eos_ptr(handle,rq,ierr)
        if (ierr /= 0) return
        call handle_allocated% assert(rq% in_use)
        if (present(gamma_melt_pt)) rq% Gamma_melt = gamma_melt_pt
        if (present(rsi_melt_pt)) rq% rsi_melt = rsi_melt_pt
        if (present(nuclide_abundance_threshold)) rq% Ythresh = nuclide_abundance_threshold
        if (present(pasta_transition_in_fm3)) rq% pasta_transition = pasta_transition_in_fm3
        if (present(cluster_transition_in_fm3)) rq% cluster_transition = cluster_transition_in_fm3
        if (present(suppress_warnings)) rq% suppress_warnings = suppress_warnings
    end subroutine dStar_eos_set_controls
    
    ! helper functions for the nuclear eos
    function nuclear_volume_fraction(rho,ionic,nuclear_radius) result(chi)
        use nucchem_def, only: composition_info_type
        use constants_def
        real(dp), intent(in) :: rho ! gram/cm**3
        type(composition_info_type), intent(in) :: ionic
        real(dp), intent(in) :: nuclear_radius  ! fm
        real(dp) :: chi
        chi = onethird*fourpi*(nuclear_radius*fm_to_cm)**3  &
            & *(rho*(1.0-ionic%Yn)/amu)
    end function nuclear_volume_fraction
    
    function neutron_wavenumber(rho,ionic,chi) result(k)
        use nucchem_def, only: composition_info_type
        use constants_def
        real(dp), intent(in) :: rho ! g/cm**3
        type(composition_info_type), intent(in) :: ionic
        real(dp), intent(in) :: chi
        real(dp) :: k   ! fm
        real(dp) :: n   ! cm**-3
        n = rho*ionic%Yn/amu/(1.0-chi)
        k = (threepisquare*n)**onethird / cm_to_fm
    end function neutron_wavenumber
    
    subroutine eval_crust_eos( &
        &   dStar_eos_handle,rho,T,ionic,ncharged,charged_ids,Yion,Tcs, &
        &   res,phase,chi,components)
        use exceptions_lib
        use nucchem_def, only: composition_info_type
        use superfluid_def, only: max_number_sf_types, neutron_1S0
        use dStar_eos_private_def
        use electron_eos
        use ion_eos
        use neutron_eos
        use radiation_eos
        use constants_def
        
        integer, intent(in) :: dStar_eos_handle
        real(dp), intent(in) :: rho,T
        type(composition_info_type), intent(in) :: ionic
        integer, intent(in) :: ncharged
        integer, dimension(ncharged), intent(in) :: charged_ids
            ! ids of the charged species
        real(dp), dimension(ncharged), intent(in) :: Yion
            ! renormalized abunances of charged species Yion = Y/(1-Yn)
        real(dp), dimension(max_number_sf_types), intent(in) :: Tcs
        real(dp), dimension(num_dStar_eos_results) :: res
        integer, intent(out) :: phase
        real(dp), intent(inout) :: chi
            ! volume fraction of nucleus; if input with value 
            ! use_default_nuclear_size, chi is computed and the new value is 
            ! returned in res.  Otherwise, the value that is input is used by 
            ! the code and unaltered.
        type(crust_eos_component), intent(out), dimension(num_crust_eos_components), optional :: components

        type(dStar_eos_general_info), pointer :: rq     
        real(dp) :: ne,rs,Gamma_e,nek, nekT,f_e, u_e, p_e, s_e, cv_e, dpr_e, dpt_e, eta_e,mu_e
        real(dp) :: f_ex, u_ex, p_ex, s_ex, cv_ex, dpr_ex, dpt_ex, uexfac, pexfac, sexfac
        real(dp) :: n_i,nik,nikT,f_i,u_i,p_i,s_i,cv_i,dpr_i,dpt_i,uifac,pifac,sifac
        real(dp) :: nn,f_n,u_n,p_n,s_n,cv_n,dpr_n,dpt_n,mu_n,unfac,pnfac,snfac,Tns
        real(dp) :: f_r,u_r,p_r,s_r,cv_r,dpr_r,dpt_r
        real(dp) :: Gamma,ionQ,p,u,s,cv,dpr,dpt,gamma3m1,gamma1,grad_ad,cp
        integer :: ierr
        type(assertion) :: neutron_volume=assertion(scope='eval_crust_eos', &
        &   message='nucleus filling fraction < 1')
        type(assertion) :: got_pointer=assertion(scope='eval_crust_eos', &
        &   message='got pointer')
        type(warning) :: ion_eos_warning=warning(scope='eval_crust_eos', &
        &   message='problem in ion eos')
        
        call dStar_eos_ptr(dStar_eos_handle, rq, ierr)
        call got_pointer% assert(ierr == 0)
        
        ! electrons...
        ! some ion quantities are defined in terms of these as well.
        ne = rho*ionic%Ye/amu
        rs = (3.0/fourpi/ne)**onethird / a_Bohr
        Gamma_e = 2.0*Rydberg/boltzmann/T/rs
        
        Gamma = 0.0
        ionQ = 0.0
                
        ! set the nuclear size
        if (chi == use_default_nuclear_size) then
            chi = nuclear_volume_fraction(rho,ionic,default_nuclear_radius)
        end if
        
        ! electrons
        call get_helm_eos_results(rho,T,ionic,f_e,u_e,p_e,s_e,cv_e,dpr_e,dpt_e,eta_e)
        call ee_exchange(Gamma_e, rs, f_ex, u_ex, p_ex, s_ex, cv_ex, dpr_ex, dpt_ex)
        nek = ne*boltzmann
        nekT = nek*T
        uexfac = nekT/rho
        pexfac = nekT
        sexfac = nek/rho
        ! u = u + u_e + u_ex*nekT/rho
        ! p = p + p_e + p_ex*nekT
        ! s = s + s_e + s_ex*nek/rho
        ! cv = cv + cv_e + cv_ex*nek/rho
        ! dpr = dpr +dpr_e + dpr_ex*nekT
        ! dpt = dpt +dpt_e + dpt_ex*nekT
!       u_ex = 0.0; p_ex = 0.0; s_ex = 0.0; cv_ex = 0.0; dpr_ex = 0.0; dpt_ex = 0.0
        
        ! ions
        call ion_mixture(rq,rs,Gamma_e,ionic,ncharged,charged_ids,Yion, &
            & Gamma,ionQ,phase,f_i,u_i,p_i,s_i,cv_i,dpr_i,dpt_i,ierr)
        
        ! ignore warning on strong quantum effects
        if (ierr == 0 .or. ierr == strong_quantum_effects) then
            n_i = ne/ionic%Z
            nik = n_i*boltzmann
            nikT = nik*T
            uifac = nikT/rho
            pifac = nikT
            sifac = nik/rho
        else if (ion_eos_warning% raised(ierr))
            uifac = 0.0; pifac = 0.0; sifac = 0.0
        end if
        
        ! local density of neutrons
        call neutron_volume% assert(chi < 1.0)
        nn = rho*ionic%Yn/amu/(1.0-chi)
        Tns = Tcs(neutron_1S0)
        call MB77(nn,T,Tns,f_n,u_n,p_n,s_n,cv_n,dpr_n,dpt_n)
        unfac = avogadro*ionic% Yn
        pnfac = 1.0
        snfac = unfac
        
        ! radiation
        call get_radiation_eos(rho,T,f_r,u_r,p_r,s_r,cv_r,dpr_r,dpt_r)

        ! stuff results into output structure
        p = p_e + p_ex*pexfac + p_i*pifac + p_n*pnfac + p_r
        u = u_e + u_ex*uexfac + u_i*uifac + u_n*unfac + u_r
        s = s_e + s_ex*sexfac + s_i*sifac + s_n*snfac + s_r
        cv = cv_e + cv_ex*sexfac + cv_i*sifac + cv_n*snfac + cv_r
        dpr = dpr_e + dpr_ex*pexfac + dpr_i*pifac + dpr_n*pnfac + dpr_r
        dpt = dpt_e + dpt_ex*pexfac + dpt_i*pifac + dpt_n*pnfac + dpt_r
        gamma3m1 = dpt/rho/T/cv
        gamma1 = (dpr + dpt*gamma3m1)/p
        grad_ad = gamma3m1/gamma1
        cp = gamma1*p*cv/dpr
        mu_e = (eta_e * boltzmann*T)*ergs_to_mev
        if (ionic% Yn < rq% Ythresh) then
            mu_n = 0.0
        else
            mu_n = (u_n -T*s_n + p_n/nn)*ergs_to_mev
        end if
        
        res(i_lnP) = log(p)
        res(i_lnE) = log(u)
        res(i_lnS)  = log(s)
        res(i_grad_ad) = grad_ad
        res(i_chiRho) = dpr/p
        res(i_chiT) = dpt/p
        res(i_Cp) = cp
        res(i_Cv) = cv
        res(i_Chi) = chi
        res(i_Gamma) = Gamma
        res(i_Theta) = 1.0/ionQ
        res(i_mu_e) = mu_e
        res(i_mu_n) = mu_n
        res(i_Gamma1) = gamma1
        res(i_Gamma3) = gamma3m1
        
        if (present(components)) then
            components(icrust_eos_ep)% P = p_e + p_ex*pexfac
            components(icrust_eos_ep)% E = u_e + u_ex*uexfac
            components(icrust_eos_ep)% S = s_e + s_ex*sexfac
            components(icrust_eos_ep)% F = f_e + f_ex*uexfac
            components(icrust_eos_ep)% Cv = cv_e + cv_ex*sexfac
            components(icrust_eos_ep)% dPdlnRho = (dpr_e + dpr_ex*pexfac)/(p_e+p_ex*pexfac)
            components(icrust_eos_ep)% dPdlnT = (dpt_e + dpt_ex*pexfac)/(p_e+p_ex*pexfac)

            components(icrust_eos_ion)% P = p_i*pifac
            components(icrust_eos_ion)% E = u_i*uifac
            components(icrust_eos_ion)% S = s_i*sifac
            components(icrust_eos_ion)% F = f_i*uifac
            components(icrust_eos_ion)% Cv = cv_i*sifac
            components(icrust_eos_ion)% dPdlnRho = dpr_i/p_i
            components(icrust_eos_ion)% dPdlnT = dpt_i/p_i

            if (ionic% Yn < rq% Ythresh) then
                components(icrust_eos_neutron) = crust_eos_component(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
            else
                components(icrust_eos_neutron)% P = p_n*pnfac
                components(icrust_eos_neutron)% E = u_n*unfac
                components(icrust_eos_neutron)% S = s_n*snfac
                components(icrust_eos_neutron)% F = f_n*unfac
                components(icrust_eos_neutron)% Cv = cv_n*snfac
                components(icrust_eos_neutron)% dPdlnRho = dpr_n/p_n
                components(icrust_eos_neutron)% dPdlnT = dpt_n/p_n
            end if
            
            components(icrust_eos_radiation)% P = p_r
            components(icrust_eos_radiation)% E = u_r
            components(icrust_eos_radiation)% S = s_r
            components(icrust_eos_radiation)% F = f_r
            components(icrust_eos_radiation)% Cv = cv_r
            components(icrust_eos_radiation)% dPdlnRho = dpr_r/p_r
            components(icrust_eos_radiation)% dPdlnT = dpt_r/p_r
        end if
    end subroutine eval_crust_eos

end module dStar_eos_lib
