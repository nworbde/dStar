module conductivity_lib
    use conductivity_def

contains

    subroutine conductivity_startup(ierr)
        use exceptions_lib
        use constants_def, only: dstar_data_dir
        use PPP_electron
        implicit none
        integer, intent(out) :: ierr
        character(len=256) :: cond_datadir
        type(alert) :: status = alert(scope='conductivity_startup')
        
        call conductivity_def_init
        cond_datadir = trim(dstar_data_dir) // '/conductivity'
        call load_PPP_electron_table(cond_datadir,ierr)
        if (ierr == unable_to_load_table) then
            call status% report( &
            &   'low-density electron conduction will be computed using analytical fmla.')
        end if
        call construct_interpolation_coefficients(ierr)
        if (ierr == unable_to_compute_interpolation) then
            call status% report( &
            &   'low-density electron conduction will be computed using analytical fmla.')
        end if
        conductivity_is_initialized = .TRUE.
    end subroutine conductivity_startup
    
    subroutine conductivity_shutdown()
        use PPP_electron
        call free_PPP_electron_table
        conductivity_is_initialized = .FALSE.
    end subroutine conductivity_shutdown
    
    function alloc_conductivity_handle(ierr) result(handle)
        use exceptions_lib
        integer, intent(out) :: ierr
        integer :: handle
        type(alert) :: status=alert(scope='alloc_conductivity_handle')
        handle = do_alloc_conductivity(ierr)
        if (ierr /= 0) then
            select case(ierr)
            case(-1)
                call status% report('unable to allocate conductivity: no free handles')
            case(-2)
                call status% report('unable to allocate conductivity: bad handle')
            end select
        end if
    end function alloc_conductivity_handle
    
    subroutine free_conductivity_handle(handle)
        integer, intent(in) :: handle
        call do_free_conductivity(handle)
    end subroutine free_conductivity_handle
    
    subroutine conductivity_set_controls(handle, include_neutrons, &
    &   include_superfluid_phonons, which_ee_scattering, which_eQ_scattering, &
    &   lgrho_rad_off, lgrho_rad_on, &
    &   max_lgrho_table, min_lgrho_table)
        use exceptions_lib
        implicit none
        integer, intent(in) :: handle
        logical, intent(in), optional :: include_neutrons
        logical, intent(in), optional :: include_superfluid_phonons
        integer, intent(in), optional :: which_ee_scattering
        integer, intent(in), optional :: which_eQ_scattering
        real(dp), intent(in), optional :: lgrho_rad_off, lgrho_rad_on
        real(dp), intent(in), optional :: max_lgrho_table, min_lgrho_table
        type(conductivity_general_info), pointer :: rq
        integer :: ierr
        type(assertion) :: got_pointer=assertion(scope='conductivity_set_controls', &
        &   message='got conductivity pointer')
        type(assertion) :: got_handle=assertion(scope='conductivity_set_controls', &
        &   message='conductivity handle is in use')
        type(alert) :: status=alert(scope='conductivity_set_controls')
        character(len=128) :: msg
        
        call get_conductivity_ptr(handle,rq,ierr)
        call got_pointer% assert(ierr == 0)
        call got_handle% assert(rq% in_use)
        
        if (present(include_neutrons)) rq% include_neutrons = include_neutrons
        if (present(include_superfluid_phonons)) rq% include_superfluid_phonons = &
        &   include_superfluid_phonons
        if (present(which_ee_scattering)) then
            if (which_ee_scattering==icond_sy06 .or. which_ee_scattering==icond_pcy) then
                rq% ee_scattering_fmla = which_ee_scattering
            else
                write(msg,'(a,i0,a)') 'unknown flag ',which_ee_scattering,' for ee scattering'
                call status% report(msg)
            end if
        end if
        if (present(which_eQ_scattering)) then
            if (which_eQ_scattering==icond_eQ_potekhin .or. which_eQ_scattering==icond_eQ_page) then
                rq% eQ_scattering_fmla = which_eQ_scattering
            else
                write(msg,'(a,i0,a)') 'unknown flag ',which_eQ_scattering,' for eQ scattering'
                call status% report(msg)
            end if
        end if
        if (present(lgrho_rad_off)) rq% rad_full_off_lgrho = lgrho_rad_off
        if (present(lgrho_rad_on)) rq% rad_full_on_lgrho = lgrho_rad_on
        if (present(max_lgrho_table)) rq% tab_off_lgrho = max_lgrho_table
        if (present(min_lgrho_table)) rq% tab_on_lgrho = min_lgrho_table
    end subroutine conductivity_set_controls

    subroutine get_thermal_conductivity(handle, &
    &   rho,T,chi,Gamma,eta,mu_e,ionic,Tcn,K)
        use exceptions_lib
        use nucchem_def, only: composition_info_type
        use eval_conductivity
        integer, intent(in) :: handle
        real(dp), intent(in) :: rho,T,Gamma,eta, mu_e, chi
        type(composition_info_type), intent(in) :: ionic
        real(dp), intent(in) :: Tcn ! neutron critical temperature
        type(conductivity_components), intent(out) :: K
        type(conductivity_general_info), pointer :: rq
        integer :: ierr
        type(assertion) :: got_pointer=assertion(scope='get_thermal_conductivity')

        call get_conductivity_ptr(handle,rq,ierr)
        call got_pointer% assert(ierr == 0)
        call conductivity(rq,rho,T,chi,Gamma,eta,mu_e,ionic,Tcn,K)
    end subroutine get_thermal_conductivity

    subroutine get_core_thermal_conductivity(nn,np,mneff,mpeff,T,Tcs,K)
        ! nn  := neutron density/fm**-3
        ! np  := proton density/fm**-3
        ! mneff := neutron effective mass/neutron rest mass
        ! mpeff := proton effective mass/proton rest mass
        ! T   := temperature (K)
        ! Tcs := critical temperatures (K)
        ! K   := thermal conductivity (cgs)
        use neutron_conductivity, only: core_neutron_conductivity
        use superfluid_def, only: max_number_sf_types
        real(dp), intent(in) :: nn,np,mneff,mpeff,T,Tcs(max_number_sf_types)
        real(dp), intent(out) :: K

        K = core_neutron_conductivity(nn,np,mneff,mpeff,T,Tcs)
    end subroutine get_core_thermal_conductivity

end module conductivity_lib
