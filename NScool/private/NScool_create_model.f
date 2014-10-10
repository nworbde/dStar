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

    integer, parameter :: number_table_pts = 128
    real(dp), parameter :: lgT_tab_min = 7.0
    real(dp), parameter :: lgT_tab_max = 10.0
    
contains
    
    subroutine do_create_crust_model(id, ierr)
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        type(NScool_info), pointer :: s
        
        call get_NScool_info_ptr(id,s,ierr)
        call do_setup_crust_zones(s, ierr)
        if (failed('do_setup_crust_zones')) return
        
        call do_setup_crust_composition(s, ierr)
        if (failed('do_setup_crust_composition')) return
        
        call do_setup_crust_transport(s, ierr)
        if (failed('do_setup_crust_transport')) return
        
    contains
        function failed(msg)
            use, intrinsic :: iso_fortran_env, only: error_unit
            character(len=*), intent(in) :: msg
            logical :: failed
            failed = .FALSE.
            if (ierr == 0) return
            write (*,*) 'do_create_crust_model failed for ',trim(msg)
            failed = .TRUE.
        end function failed
    end subroutine do_create_crust_model
    
    subroutine do_setup_crust_zones(s, ierr)
        use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
        use constants_def
        use nucchem_lib
        use superfluid_lib
        use dStar_eos_lib
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
        
        write (error_unit,*) 'establishing crust zones...'
        
        write (error_unit,*) 'initializing microphysics...'
        call nucchem_init(trim(dStar_data_dir),ierr)
        if (failure('nucchem_init')) return
	
        call sf_startup(trim(dStar_data_dir),ierr)
        if (failure('sf_startup')) return
	
        call sf_load_gaps(trim(s% which_proton_1S0_gap), trim(s% which_neutron_1S0_gap), &
            & trim(s% which_neutron_3P2_gap), ierr)
!         call sf_load_gaps('ns','sfb03','t72',ierr)
        if (failure('sf_load_gaps')) return
	
        call dStar_eos_startup(trim(dStar_data_dir))
        if (failure('dStar_eos_startup')) return
	
        s% eos_handle = alloc_dStar_eos_handle(ierr)
        if (failure('alloc_dStar_eos_handle')) return
    
        ! switch off the warnings about quantum effects
        call dStar_eos_set_controls(s% eos_handle,suppress_warnings=.TRUE.)
        if (s% eos_gamma_melt_pt > 0.0)  &
        & call dStar_eos_set_controls(s% eos_handle,gamma_melt_pt=s% eos_gamma_melt_pt)
        if (s% eos_rsi_melt_pt > 0.0)  &
        & call dStar_eos_set_controls(s% eos_handle,rsi_melt_pt=s% eos_rsi_melt_pt)
        if (s% eos_nuclide_abundance_threshold > 0.0)  &
        & call dStar_eos_set_controls(s% eos_handle,nuclide_abundance_threshold=s% eos_nuclide_abundance_threshold)
        if (s% eos_pasta_transition_in_fm3 > 0.0)  &
        & call dStar_eos_set_controls(s% eos_handle,pasta_transition_in_fm3=s% eos_pasta_transition_in_fm3)
        if (s% eos_cluster_transition_in_fm3 > 0.0)  &
        & call dStar_eos_set_controls(s% eos_handle,cluster_transition_in_fm3=s% eos_cluster_transition_in_fm3)
            
        write (error_unit,*) 'loading crust model...'
        call dStar_crust_startup(trim(dStar_data_dir),ierr)
        if (failure('dStar_crust_startup')) return
        
        call dStar_crust_load_table('hz90',s% eos_handle, s% Tcore, ierr)
        if (failure('dStar_crust_load_table')) return

    	call dStar_atm_startup(trim(dStar_data_dir),ierr)
        if (failure('dStar_atm_startup')) return

        write(error_unit,*) 'integrating TOV equations...'
        allocate(y(num_tov_variables))

        call tov_integrate(log10(s% Pcore), log10(s% Psurf), s% Mcore, s% Rcore, s% target_resolution_lnP, y, ierr)
        if (failure('tov_integrate')) return
        deallocate(y)
        
        ! allocate the info arrays
        stov => tov_model
        s% nz = stov% nzs - 1
        s% nisos = dStar_crust_get_composition_size()
        call allocate_NScool_info_arrays(s, ierr)
        
        ! facial information
        s% P_bar(1:s% nz) = stov% pressure(stov% nzs:2:-1) * pressure_g
        s% ePhi_bar(1:s% nz) = exp(stov% potential(stov% nzs:2:-1))
        s% e2Phi_bar(1:s% nz) = s% ePhi_bar(1:s% nz)**2
        s% ePhicore = exp(stov% potential(1))
        s% m(1:s% nz) = stov% baryon(stov% nzs:2:-1) * mass_g
        s% area(1:s% nz) = fourpi*(stov% radius(stov% nzs:2:-1)*length_g + s% Rcore*1.0e5)**2
        s% eLambda_bar(1:s% nz) = 1.0/sqrt(1.0-2.0*(stov% mass(stov% nzs:2:-1)+s% Mcore)/ &
        &   (stov% radius(stov% nzs:2:-1) + s% Rcore*1.0e5/length_g))
        
        ! body information
        s% dm(1:s% nz-1) = s% m(1:s% nz-1) - s% m(2:s% nz)
        s% dm(s% nz) = s% m(s% nz)
        
        ! mean density of a cell
        s% rho(1:s% nz) = s% dm(1:s% nz)/(stov% volume(stov% nzs:2:-1) - stov% volume(stov% nzs-1:1:-1))/length_g**3
        
        ! interpolate metric functions
        s% eLambda(1:s% nz-1) = 0.5*(s% eLambda_bar(1:s% nz-1) + s% eLambda_bar(2:s% nz))
        s% eLambda(s% nz) = 0.5*(s% eLambda_bar(s% nz) +  &
        &   1.0/sqrt(1.0-2.0*s% Mcore/(s% Rcore*1.0e5/length_g)))
        s% ePhi(1:s% nz-1) = 0.5*(s% ePhi_bar(1:s% nz-1) + s% ePhi_bar(2:s% nz))
        s% ePhi(s% nz) = 0.5*(s% ePhi_bar(s% nz) + s% ePhicore)
        s% e2Phi(1:s% nz) = s% ePhi(1:s% nz)**2
        
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
        
        ! now set the surface gravity and load the atmosphere
        s% grav = stov% mass(stov% nzs) + s% Mcore * mass_g * Gnewton /  &
        &   (stov% radius(stov% nzs)*length_g + s% Rcore*1.0e5)**2 * s% eLambda_bar(1)

        Plight = s% grav * 10.0_dp**s% lg_atm_light_element_column
    	call dStar_atm_load_table('pcy97',s% grav, Plight,ierr)
        if (failure('dStar_atm_load_table')) return

    contains
        function failure(str)
            character(len=*), intent(in) :: str
            logical :: failure
            failure = .FALSE.
            if (ierr == 0) return
            write (error_unit,*) 'do_setup_crust_zones failure in ',trim(str)
            failure = .TRUE.
        end function failure
        
    end subroutine do_setup_crust_zones
    
    subroutine do_setup_crust_composition(s, ierr)
        use, intrinsic :: iso_fortran_env, only: error_unit
        use nucchem_def
        use nucchem_lib
        use storage

        type(NScool_info), pointer :: s
        integer, intent(out) :: ierr
        real(dp), allocatable, dimension(:) :: lgP_bar
        
        ! this routine assumes that we have already run do_setup_crust_zones
        ierr = 0
        allocate(lgP_bar(s% nz))
        lgP_bar = log10(s% P_bar)
        
        call allocate_NScool_iso_arrays(s, ierr)
        call dStar_crust_get_composition(lgP_bar, s% ncharged, s% charged_ids, s% Yion_bar, s% Xneut_bar, s% ionic_bar, ierr)
        if (ierr /= 0) return
        
        ! can fix the impurity parameter to a specified value; this should be replaced by
        ! a function call to allow a more customizable option
        if (s% fix_Qimp) then
            write (error_unit,'(a,f7.2)') 'setting Qimp = ',s% Qimp
            s% ionic_bar(1:s% nz)% Q = s% Qimp
        end if

        ! for the cells, *for now*, inherit composition of the top face
        s% Yion(1:s% ncharged, 1:s% nz) = s% Yion_bar(1:s% ncharged, 1:s% nz)
        s% ionic(1:s% nz) = s% ionic_bar(1:s% nz)
        s% Xneut(1:s% nz) = s% Xneut_bar(1:s% nz)
        
        deallocate(lgP_bar)
    end subroutine do_setup_crust_composition
    
    subroutine do_setup_crust_transport(s, ierr)
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
        use utils_lib, only: is_bad_num
        
        type(NScool_info), pointer :: s
        integer, intent(out) :: ierr
        logical, dimension(num_crust_nu_channels) :: nu_channels
        logical, dimension(num_conductivity_channels) :: cond_channels
        type(crust_neutrino_emissivity_channels) :: eps_nu
        integer :: iz, itemp, ieos
        integer :: eos_phase
        real(dp), dimension(num_dStar_eos_results) :: eos_results
        type(conductivity_components) :: Kcomponents
        real(dp) :: chi, Ttab, enu_tab
        real(dp), dimension(:), pointer :: work=>null()
        real(dp), dimension(:,:), pointer :: lnEnu_val, lnKcond_val, lnCp_val
        real(dp), dimension(:), pointer :: lnEnu_interp, lnKcond_interp, lnCp_interp
        real(dp) :: nn, kn, Tc(max_number_sf_types)
        type(crust_eos_component), dimension(num_crust_eos_components) :: components
        ! for error checking
        real(dp), dimension(:), allocatable :: delP
        
        ierr = 0
        s% n_tab = number_table_pts
        
        call allocate_NScool_work_arrays(s, ierr)
        if (ierr /= 0) return
        
        ! set up the conductivity and neutrino channels
        nu_channels = [ s% use_crust_nu_pair, s% use_crust_nu_photo, s% use_crust_nu_plasma, &
        & s% use_crust_nu_bremsstrahlung, s% use_crust_nu_pbf ]
        
        cond_channels = [ s% use_ee_conductivity, &
        & s% use_ei_conductivity, s% use_eQ_conductivity, s% use_sf_conductivity ]
        
        s% tab_lnT(1:s% n_tab) = [(lgT_tab_min*ln10 + (lgT_tab_max-lgT_tab_min)*ln10*real(itemp-1,dp)/real(s% n_tab-1,dp), &
        &   itemp = 1, s% n_tab)]
        
        ! cell average quantities: Cp and enu
        do iz = 1, s% nz
            lnCp_val(1:4,1:s% n_tab) => s% tab_lnCp(1:4*s% n_tab, iz)
            lnEnu_val(1:4,1:s% n_tab) => s% tab_lnEnu(1:4*s% n_tab, iz)
            do itemp = 1, s% n_tab
                Ttab = exp(s% tab_lnT(itemp))
                chi = use_default_nuclear_size
                call eval_crust_eos( &
                &   s% eos_handle, s% rho(iz), Ttab, s% ionic(iz), s% ncharged, s% charged_ids, s% Yion(1:s% ncharged,iz), &
                &   eos_results, eos_phase, chi, components)
                lnCp_val(1,itemp) = log(eos_results(i_Cp))
                if (is_bad_num(lnCp_val(1,itemp)) .and. itemp == 40) then
                    print *,'bad CV'
                    do ieos = 1, num_crust_eos_components
                        print *, ieos
                        print *,components(ieos)
                    end do
                    print *, s% ionic(iz), s% Yion(1:s% ncharged,iz)
                end if
                if (itemp == 1) s% P(iz) = exp(eos_results(i_lnP))
                
                ! perform interpolation of density to first cell face
                if (iz == 1 .and. itemp == 1) then
                    s% rho_bar(iz) = s% rho(iz)*(s% P_bar(iz)/s% P(iz))**(1.0/eos_results(i_chiRho))
                end if
                
                nn = s% rho(iz)*s% ionic(iz)% Yn/amu/(1.0-chi)
                kn = (1.5*pi**2*nn)**onethird / cm_to_fm
                call sf_get_results(0.0_dp,kn,Tc)
                
                call get_crust_neutrino_emissivity(s% rho(iz), Ttab, s% ionic(iz), chi, Tc(neutron_1S0),  &
                &   eps_nu, nu_channels)
                lnEnu_val(1,itemp) = log(eps_nu% total/s% rho(iz))
                
            end do

            ! add in the shell Urca cooling
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
            lnCp_interp(1:4*s% n_tab) => s% tab_lnCp(1:4*s% n_tab,iz)
            call interp_pm(s% tab_lnT, s% n_tab, lnCp_interp, pm_work_size, work, &
            &   'do_setup_crust_transport: lnCp', ierr)            
            deallocate(work)
            
        end do
        
        ! facial quantities: Kcond
        ! compute dp for error checking
        allocate(delP(s% nz))
        do iz = 1, s% nz
            lnKcond_val(1:4,1:s% n_tab) => s% tab_lnK(1:4*s% n_tab, iz)
            chi = use_default_nuclear_size
            do itemp = 1, s% n_tab
                Ttab = exp(s% tab_lnT(itemp))
                call eval_crust_eos( &
                &   s% eos_handle, s% rho_bar(iz), Ttab, s% ionic_bar(iz), s% ncharged, s% charged_ids,  &
                &   s% Yion_bar(1:s% ncharged,iz), eos_results, eos_phase, chi)
                if (itemp == 1) delP(iz) = abs(1.0 - exp(eos_results(i_lnP))/s% P_bar(iz))
                                
                call get_thermal_conductivity(s% rho_bar(iz), Ttab, chi, eos_results(i_Gamma),  &
                &   eos_results(i_Theta), s% ionic_bar(iz), &
                &   Kcomponents, use_pcy=s% use_pcy_for_ee_scattering, use_page=s% use_page_for_eQ_scattering, &
                &   which_components=cond_channels)
                lnKcond_val(1,itemp) = log(Kcomponents% total)

            end do
            allocate(work(s% n_tab*pm_work_size))
            lnKcond_interp(1:4*s% n_tab) => s% tab_lnK(1:4*s% n_tab,iz)
            call interp_pm(s% tab_lnT, s% n_tab, lnKcond_interp, pm_work_size, work, &
        &   'do_setup_crust_transport: lnK', ierr)
            deallocate(work)
        end do
!         write (*,'(a,es15.8,a,es15.8)') 'max dP = ',maxval(delP),' at ',s% P_bar(maxloc(delP))
        deallocate(delP)
        
    end subroutine do_setup_crust_transport

end module create_model
