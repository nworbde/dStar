module create_model
     use NScool_def
    use constants_def, only: dp
    use constants_lib
    use nucchem_def
    use nucchem_lib
    use superfluid_lib
    use dStar_eos_lib
    use dStar_crust_def
    use dStar_crust_lib
    use dStar_atm_def
    use dStar_atm_lib
    use NScool_crust_tov

    ! for tabulation of coefficients as a fcn of temprerature
    integer, parameter :: number_table_pts = 128
    real(dp), parameter :: lgT_tab_min = 7.0
    real(dp), parameter :: lgT_tab_max = 10.0
    
    character(len=*), parameter, private :: indent1 = '(tr4,"* ",a)'
    
contains
    
    subroutine do_create_crust_model(id, ierr)
        use exceptions_lib
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        type(NScool_info), pointer :: s
        type(assertion) :: got_pointer=assertion(scope='do_create_crust_model')
        type(failure) :: microphysics_error=failure(scope='do_create_crust_model', &
        &   message='staring microphysics')
        type(failure) :: zone_error=failure(scope='do_create_crust_model', &
        &   message='setting crust zones')
        type(failure) :: composition_error= &
        &   failure(scope='do_create_crust_model', &
        &   message='setting crust composition')
        type(failure) :: transport_error=failure( &
        &   scope='do_create_crust_model',message='setting crust transport')
        
        call get_NScool_info_ptr(id,s,ierr)
        call got_pointer% assert(ierr==0)
        
        call do_startup_microphysics(s, ierr)
        if (microphysics_error% raised(ierr)) return
        
        call do_setup_crust_zones(s, ierr)
        if (zone_error% raised(ierr)) return
        
        call do_setup_crust_composition(s, ierr)
        if (composition_error% raised(ierr)) return
        
        call do_setup_crust_transport(s, ierr)
        if (transport_error% raised(ierr)) return
    end subroutine do_create_crust_model
    
    subroutine do_startup_microphysics(s, ierr)
        use exceptions_lib
        use constants_def
        use nucchem_lib
        use superfluid_lib
        use dStar_eos_lib
        use conductivity_def
        use conductivity_lib
        use NScool_private_def, only: dStar_data_dir

        type(NScool_info), pointer :: s
        integer, intent(out) :: ierr
        type(failure) :: nucchem_error=failure(scope='do_startup_microphysics', &
        &   message='initializing nucchem module')
        type(failure) :: superfluid_error=failure(scope='do_startup_microphysics', &
        &   message='starting superfluid module')
        type(failure) :: gap_error=failure(scope='do_startup_microphysics', &
        &   message='loading gaps')
        type(failure) :: eos_error=failure(scope='do_startup_microphysics', &
        &   message='starting eos module')
        type(failure) :: cond_error=failure(scope='do_startup_microphysics', &
        &   message='starting conductivity module')
        type(alert) :: status=alert(scope='do_startup_microphysics')
        
        call status% report('Initializing nuclides')
        call nucchem_init(trim(dStar_data_dir),ierr)
        if (nucchem_error% raised(ierr)) return
        
        call status% report('Loading superfluid gaps')
        call sf_startup(trim(dStar_data_dir),ierr)
        if (superfluid_error% raised(ierr)) return
        call sf_load_gaps( &
        &   trim(s% which_proton_1S0_gap),  &
        &   trim(s% which_neutron_1S0_gap), &
        &   trim(s% which_neutron_3P2_gap), ierr)
        if (gap_error% raised(ierr)) return
        sf_scale(1:max_number_sf_types) = s% scale_sf_critical_temperatures
    
        call status% report('Setting EOS options')
        call dStar_eos_startup(trim(dStar_data_dir))
        s% eos_handle = alloc_dStar_eos_handle(ierr)
        if (eos_error% raised(ierr)) return
        ! switch off the warnings about quantum effects
        call dStar_eos_set_controls(s% eos_handle,suppress_warnings=.TRUE.)
        if (s% eos_gamma_melt_pt > 0.0)  &
        &   call dStar_eos_set_controls(s% eos_handle, &
        &   gamma_melt_pt=s% eos_gamma_melt_pt)
        if (s% eos_rsi_melt_pt > 0.0)  &
        &   call dStar_eos_set_controls(s% eos_handle, &
        &   rsi_melt_pt=s% eos_rsi_melt_pt)
        if (s% eos_nuclide_abundance_threshold > 0.0)  &
        & call dStar_eos_set_controls(s% eos_handle, &
        &   nuclide_abundance_threshold=s% eos_nuclide_abundance_threshold)
        if (s% eos_pasta_transition_in_fm3 > 0.0)  &
        & call dStar_eos_set_controls(s% eos_handle, &
        &   pasta_transition_in_fm3=s% eos_pasta_transition_in_fm3)
        if (s% eos_cluster_transition_in_fm3 > 0.0)  &
        & call dStar_eos_set_controls(s% eos_handle, &
        &   cluster_transition_in_fm3=s% eos_cluster_transition_in_fm3)
 
        call status% report('Setting thermal conductivity options')
        call conductivity_startup(trim(dStar_data_dir))
        s% cond_handle = alloc_conductivity_handle(ierr)
        if (cond_error% raised(ierr)) return
        call conductivity_set_controls(s% cond_handle, &
        &   include_neutrons=s% use_neutron_conductivity, &
        &   include_superfluid_phonons=s% use_superfluid_phonon_conductivity)
        if (s% rad_full_off_lgrho > 0.0) &
        &   call conductivity_set_controls(s% cond_handle, &
        &   lgrho_rad_off=s% rad_full_off_lgrho)
        if (s% rad_full_on_lgrho > 0.0) &
        &   call conductivity_set_controls(s% cond_handle, &
        &   lgrho_rad_on=s% rad_full_on_lgrho)
        if (s% use_pcy_for_ee_scattering)  &
        &   call conductivity_set_controls(s% cond_handle, &
        &   which_ee_scattering=icond_pcy)
        if (s% use_page_for_eQ_scattering)  &
        &   call conductivity_set_controls(s% cond_handle, &
        &   which_eQ_scattering=icond_eQ_page)            
    end subroutine do_startup_microphysics
    
    subroutine do_setup_crust_zones(s, ierr)
        use exceptions_lib
        use constants_def
        use nucchem_lib
        use superfluid_lib
        use dStar_eos_lib
        use conductivity_def
        use conductivity_lib
        use dStar_crust_def
        use dStar_crust_lib
        use NScool_private_def, only: dStar_data_dir
        use NScool_crust_tov
        use storage
        
        type(NScool_info), pointer :: s
        type(tov_model_type), pointer :: stov
        integer, intent(out) :: ierr
        real(dp), dimension(:), pointer :: y
        real(dp) :: Plight
        type(alert) :: status=alert(scope='do_setup_crust_zones')
        type(failure) :: crust_error=failure(scope='do_setup_crust_zones', &
        &   message='starting crust module')
        type(failure) :: crust_table_error=failure(scope='do_setup_crust_zones', &
        &   message='setting crust zones')
        type(failure) :: atm_error=failure(scope='do_setup_crust_zones', &
        &   message='starting atmosphere module')
        type(failure) :: tov_error=failure(scope='do_setup_crust_zones', &
        &   message='integrating TOV eqns')
        type(failure) :: atm_load_error=failure(scope='do_setup_crust_zones', &
        &   message='loading atm table')
        type(assertion) :: allocation_okay=assertion(scope='do_setup_crust_zones', &
        &   message='memory allocated')
        
        call status% report('Initializing crust')
        call dStar_crust_startup(trim(dStar_data_dir),ierr)
        if (crust_error% raised(ierr)) return

        call status% report('Initializing atmosphere')
        call dStar_atm_startup(trim(dStar_data_dir),ierr)
        if (atm_error% raised(ierr)) return

        call status% report('loading crust table')
        call dStar_crust_load_table(s% crust_composition,s% eos_handle, s% crust_reference_temperature, ierr)
        if (crust_error% raised(ierr)) return

        call status% report('Integrating TOV equations')
        allocate(y(num_tov_variables))
        call tov_integrate(log10(s% Pcore), log10(s% Psurf), s% Mcore, s% Rcore, s% target_resolution_lnP, y, ierr)
        if (tov_error% raised(ierr)) return
        deallocate(y)        
        stov => tov_model
        s% nz = stov% nzs - 1
        s% nisos = dStar_crust_get_composition_size()
        call allocate_NScool_info_arrays(s, ierr)
        call allocation_okay% assert(ierr==0)
        
        call status% report('Computing facial quantities')
        s% P_bar(1:s% nz) = stov% pressure(stov% nzs:2:-1) * pressure_g
        s% ePhi_bar(1:s% nz) = exp(stov% potential(stov% nzs:2:-1))
        s% e2Phi_bar(1:s% nz) = s% ePhi_bar(1:s% nz)**2
        s% ePhicore = exp(stov% potential(1))
        s% m(1:s% nz) = stov% baryon(stov% nzs:2:-1) * mass_g
        s% area(1:s% nz) = fourpi*(stov% radius(stov% nzs:2:-1)*length_g + s% Rcore*1.0e5)**2
        s% eLambda_bar(1:s% nz) = 1.0/sqrt(1.0-2.0*(stov% mass(stov% nzs:2:-1)+s% Mcore)/ &
        &   (stov% radius(stov% nzs:2:-1) + s% Rcore*1.0e5/length_g))
        
        call status% report('Computing zonal quantities')
        s% dm(1:s% nz-1) = s% m(1:s% nz-1) - s% m(2:s% nz)
        s% dm(s% nz) = s% m(s% nz)
        ! mean density of a cell
        s% rho(1:s% nz) = s% dm(1:s% nz)/(stov% volume(stov% nzs:2:-1) - stov% volume(stov% nzs-1:1:-1))/length_g**3
        
        call status% report('Interpolating metric functions')
        s% eLambda(1:s% nz-1) = 0.5*(s% eLambda_bar(1:s% nz-1) + s% eLambda_bar(2:s% nz))
        s% eLambda(s% nz) = 0.5*(s% eLambda_bar(s% nz) +  &
        &   1.0/sqrt(1.0-2.0*s% Mcore/(s% Rcore*1.0e5/length_g)))
        s% ePhi(1:s% nz-1) = 0.5*(s% ePhi_bar(1:s% nz-1) + s% ePhi_bar(2:s% nz))
        s% ePhi(s% nz) = 0.5*(s% ePhi_bar(s% nz) + s% ePhicore)
        s% e2Phi(1:s% nz) = s% ePhi(1:s% nz)**2
        s% P(1:s% nz-1) = 0.5*(s% P_bar(1:s% nz-1)+s% P_bar(2:s% nz))
        
        ! temperatures: set to be isothermal (exp(Phi)*T = const)
        s% T(1:s% nz) = s% Tcore * s% ePhicore / s% ePhi(1:s% nz)
        s% T_bar(1:s% nz) = s% Tcore * s% ePhicore / s% ePhi_bar(1:s%nz)
        s% lnT(1:s% nz) = log(s% T(1:s% nz))
        s% lnT_bar(1:s% nz) = log(s% T_bar(1:s% nz))
         
        ! now interpolate dm and rho to the facial points
        s% dm_bar(1) = 0.5*s% dm(1)
        s% dm_bar(2:s% nz) = 0.5*(s% dm(1:s% nz-1) + s% dm(2:s% nz))

        ! this leaves rho_bar(1) undefined... we can use Taylor expansion
        s% rho_bar(2:s% nz) = 0.5*(s% rho(1:s% nz-1)*s% dm(2:s% nz) + s% rho(2:s% nz)*s% dm(1:s% nz-1)) / &
        &   s% dm_bar(2:s% nz)
        
        call status% report('Loading atmosphere')
        ! mass in Msun
        s% Mtotal = stov% mass(stov% nzs) + s% Mcore
        ! R in km
        s% Rtotal = stov% radius(stov% nzs)*length_g*1.0e-5 + s% Rcore
        s% grav = (s% Mtotal) * mass_g * Gnewton / (s% Rtotal*1.0e5)**2 * s% eLambda_bar(1)
        Plight = s% grav * 10.0_dp**s% lg_atm_light_element_column
        call dStar_atm_load_table(s% atm_model, s% grav, Plight, s% Psurf, ierr)
        if (atm_load_error% raised(ierr)) return        
    end subroutine do_setup_crust_zones
    
    subroutine do_setup_crust_composition(s, ierr)
        use exceptions_lib
        use nucchem_def
        use nucchem_lib
        use storage

        type(NScool_info), pointer :: s
        integer, intent(out) :: ierr
        real(dp), allocatable, dimension(:) :: lgP_bar
        type(alert) :: qimp=alert(scope='do_setup_crust_composition')
        character(len=128) :: qimp_msg
        type(alert) :: status=alert(scope='do_setup_crust_composition')
                
        ! this routine assumes that we have already run do_setup_crust_zones
        ierr = 0
        allocate(lgP_bar(s% nz))
        lgP_bar = log10(s% P_bar)
        
        call status% report('Setting composition')
        call allocate_NScool_iso_arrays(s, ierr)
        call dStar_crust_get_composition_info(lgP_bar, s% ncharged, s% charged_ids, s% Yion_bar, s% Xneut_bar, s% ionic_bar, ierr)
        if (ierr /= 0) return
        
        ! can fix the impurity parameter to a specified value
        if (s% fix_Qimp) then
            write (qimp_msg,'(a,f7.2)') 'Setting Qimp = ',s% Qimp
            call qimp% report(qimp_msg)
            s% ionic_bar(1:s% nz)% Q = s% Qimp
        end if
        
        ! if there is a customized option, use that
        if (s% use_other_set_Qimp) then
            call qimp% report('Using user-defined Qimp')
            call s% other_set_Qimp(s% id, ierr)
        end if

        ! for the cells, *for now*, inherit composition of the top face
        s% Yion(1:s% ncharged, 1:s% nz) = s% Yion_bar(1:s% ncharged, 1:s% nz)
        s% ionic(1:s% nz) = s% ionic_bar(1:s% nz)
        s% Xneut(1:s% nz) = s% Xneut_bar(1:s% nz)
        
        deallocate(lgP_bar)
    end subroutine do_setup_crust_composition
    
    subroutine do_setup_crust_transport(s, ierr)
        use exceptions_lib
        use constants_def
        use nucchem_lib
        use superfluid_def
        use superfluid_lib
        use dStar_eos_def
        use dStar_eos_lib
        use conductivity_def
        use conductivity_lib
        use neutrino_def
        use neutrino_lib
        use interp_1d_def
        use interp_1d_lib
        use storage
        
        type(NScool_info), pointer :: s
        integer, intent(out) :: ierr
        logical, dimension(num_crust_nu_channels) :: nu_channels
        type(crust_neutrino_emissivity_channels) :: eps_nu
        integer :: iz, itemp, ieos
        integer :: eos_phase
        real(dp), dimension(num_dStar_eos_results) :: eos_results
        type(conductivity_components) :: Kcomponents
        real(dp) :: chi, Ttab, enu_tab
        real(dp), dimension(:), pointer :: work=>null()
        real(dp), dimension(:,:), pointer :: lnEnu_val, lnKcond_val, lnCp_val, lnGamma_val
        real(dp), dimension(:), pointer :: lnEnu_interp, lnKcond_interp, lnCp_interp, lnGamma_interp
        real(dp) :: kn, kp, Tc(max_number_sf_types)
        type(crust_eos_component), dimension(num_crust_eos_components) :: components
        ! for error checking
        real(dp), dimension(:), allocatable :: delP
        type(alert) :: status=alert(scope='do_setup_crust_transport')
        type(assertion) :: interpolation_okay=assertion(scope='do_setup_crust_transport', &
        &    message='interpolation coefficients computed')
        type(assertion) :: allocation_okay=assertion(scope='do_setup_crust_transport', &
        &    message='allocation okay')
        
        ierr = 0
        s% n_tab = number_table_pts
        
        call allocate_NScool_work_arrays(s, ierr)
        call allocation_okay% assert(ierr == 0)
        
        ! set up the conductivity and neutrino channels
        nu_channels = [ s% use_crust_nu_pair,  &
        &               s% use_crust_nu_photo,  &
        &               s% use_crust_nu_plasma, &
        &               s% use_crust_nu_bremsstrahlung,  &
        &               s% use_crust_nu_pbf ]
                
        s% tab_lnT(1:s% n_tab) = [(lgT_tab_min*ln10 + (lgT_tab_max-lgT_tab_min)*ln10*real(itemp-1,dp)/real(s% n_tab-1,dp), &
        &   itemp = 1, s% n_tab)]
        
        call status% report('Computing specific heat, neutrino emissivity tables')
        do iz = 1, s% nz
            lnCp_val(1:4,1:s% n_tab) => s% tab_lnCp(1:4*s% n_tab, iz)
            lnGamma_val(1:4,1:s% n_tab) => s% tab_lnGamma(1:4*s% n_tab, iz)
            lnEnu_val(1:4,1:s% n_tab) => s% tab_lnEnu(1:4*s% n_tab, iz)
            chi = nuclear_volume_fraction(s% rho(iz),s% ionic(iz), &
            &   default_nuclear_radius)
            kn = neutron_wavenumber(s% rho(iz),s% ionic(iz),chi)
            kp = 0.0_dp
            if (.not. s% use_other_sf_critical_temperatures) then
                call sf_get_results(kp,kn,Tc)
            else
                call s% other_sf_get_results(s% id,kp,kn,Tc)
            end if
            do itemp = 1, s% n_tab
                Ttab = exp(s% tab_lnT(itemp))
                call eval_crust_eos( &
                &   s% eos_handle, s% rho(iz), Ttab, s% ionic(iz),  &
                &   s% ncharged, s% charged_ids, s% Yion(1:s% ncharged,iz), &
                &   Tc, eos_results, eos_phase, chi, components)
                lnCp_val(1,itemp) = log(eos_results(i_Cp))
                lnGamma_val(1,itemp) = log(eos_results(i_Gamma))

                ! perform interpolation of density to first cell face
                if (iz == 1 .and. itemp == 1) then
                    s% rho_bar(iz) = s% rho(iz)*(s% P_bar(iz)/s% P(iz))**(1.0/eos_results(i_chiRho))
                end if
                                
                call get_crust_neutrino_emissivity(s% rho(iz), Ttab,  &
                &   s% ionic(iz), chi, Tc(neutron_1S0),  &
                &   eps_nu, nu_channels)
                lnEnu_val(1,itemp) = log(eps_nu% total/s% rho(iz))
                
            end do

            ! add in the shell Urca cooling: this should be a call to a (potentially user-defined) function
            if (s% turn_on_shell_Urca .and. iz < s% nz) then
                if (log10(s% P_bar(iz)) < s% lgP_shell_Urca .and. log10(s% P_bar(iz+1)) > s% lgP_shell_Urca) then
                    lnEnu_val(1,1:s% n_tab) = log(exp(lnEnu_val(1,1:s% n_tab))  &
                    &   + 1.0e36_dp*s% shell_Urca_luminosity_coeff * (exp(s% tab_lnT(1:s% n_tab))/5.0e8_dp)**5/s% dm(iz))
                end if
            end if
            
            allocate(work(s% n_tab*pm_work_size))
            lnEnu_interp(1:4*s% n_tab) => s% tab_lnEnu(1:4*s% n_tab,iz)
            call interp_pm(s% tab_lnT, s% n_tab, lnEnu_interp, pm_work_size, work, &
            &   'do_setup_crust_transport: lnEnu', ierr)
            call interpolation_okay% assert(ierr == 0)
            lnCp_interp(1:4*s% n_tab) => s% tab_lnCp(1:4*s% n_tab,iz)
            call interp_pm(s% tab_lnT, s% n_tab, lnCp_interp, pm_work_size, work, &
            &   'do_setup_crust_transport: lnCp', ierr)
            call interpolation_okay% assert(ierr == 0)
            lnGamma_interp(1:4*s% n_tab) => s% tab_lnGamma(1:4*s% n_tab, iz)
            call interp_pm(s% tab_lnT, s% n_tab, lnGamma_interp, pm_work_size, work,  &
            &   'do_setup_crust_transport: Gamma', ierr)
            call interpolation_okay% assert(ierr == 0)
            deallocate(work)
            
        end do
        
        call status% report('Computing thermal conductivity tables')
        ! compute dp for error checking
        allocate(delP(s% nz))
        do iz = 1, s% nz
            lnKcond_val(1:4,1:s% n_tab) => s% tab_lnK(1:4*s% n_tab, iz)
            do itemp = 1, s% n_tab
                Ttab = exp(s% tab_lnT(itemp))
                chi = nuclear_volume_fraction(s% rho_bar(iz),s% ionic_bar(iz),&
                &   default_nuclear_radius)
                kn = neutron_wavenumber(s% rho_bar(iz), s% ionic_bar(iz),chi)
                kp = 0.0_dp
                if (.not. s% use_other_sf_critical_temperatures) then
                    call sf_get_results(kp,kn,Tc)
                else
                    call s% other_sf_get_results(s% id,kp,kn,Tc)
                end if
                call eval_crust_eos( &
                &   s% eos_handle, s% rho_bar(iz), Ttab, s% ionic_bar(iz),  &
                &   s% ncharged, s% charged_ids, &
                &   s% Yion_bar(1:s% ncharged,iz), Tc, eos_results,  &
                &   eos_phase, chi)
                if (itemp == 1)  &
                &   delP(iz) = abs(1.0 - exp(eos_results(i_lnP))/s% P_bar(iz))
                                
                call get_thermal_conductivity(s% cond_handle, &
                &   s% rho_bar(iz), Ttab, chi, &
                &   eos_results(i_Gamma), eos_results(i_Theta), &
                &   eos_results(i_mu_e), s% ionic_bar(iz), Tc(neutron_1S0), &
                &   Kcomponents)
                lnKcond_val(1,itemp) = log(Kcomponents% total)

            end do
            allocate(work(s% n_tab*pm_work_size))
            lnKcond_interp(1:4*s% n_tab) => s% tab_lnK(1:4*s% n_tab,iz)
            call interp_pm(s% tab_lnT, s% n_tab, lnKcond_interp, pm_work_size, work, &
            &   'do_setup_crust_transport: lnK', ierr)
            call interpolation_okay% assert(ierr == 0)
            deallocate(work)
        end do
        deallocate(delP)
        
    end subroutine do_setup_crust_transport

end module create_model
