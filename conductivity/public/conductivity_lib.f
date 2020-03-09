module conductivity_lib
    use conductivity_def

contains

    subroutine conductivity_startup(datadir)
        use iso_fortran_env, only: error_unit
        use PPP_electron
        implicit none
        character(len=*), intent(in) :: datadir
        character(len=256) :: cond_datadir
        integer :: ierr
        call conductivity_def_init
        
        cond_datadir = trim(datadir) // '/conductivity'
        call load_PPP_electron_table(cond_datadir,ierr)
        if (ierr == unable_to_load_table) then
            write(error_unit, *) 'low-density electron conduction will not be computed correctly'
        end if
        call construct_interpolation_coefficients(ierr)
        if (ierr == unable_to_compute_interpolation) then
            write(error_unit, *) 'low-density electron conduction will not be computed correctly'
            PPP_tbl% is_loaded = .FALSE.
        end if
    end subroutine conductivity_startup
    
    subroutine conductivity_shutdown()
        use PPP_electron
        call free_PPP_electron_table
    end subroutine conductivity_shutdown
    
    function alloc_conductivity_handle(ierr) result(handle)
        use, intrinsic :: iso_fortran_env, only: error_unit
        integer, intent(out) :: ierr
        integer :: handle
        handle = do_alloc_conductivity(ierr)
        if (ierr /= 0) then
            select case(ierr)
            case(-1)
                write(error_unit,*) 'unable to allocate conductivty: no free handles'
            case(-2)
                write(error_unit,*) 'unable to allocate conductivty: bad handle'
            end select
        end if
    end function alloc_conductivity_handle
    
    subroutine free_conductivity_handle(handle)
        integer, intent(in) :: handle
        call do_free_conductivity(handle)
    end subroutine free_conductivity_handle
    
    subroutine conductivity_set_controls(handle,include_electrons, include_neutrons, &
    &   include_superfluid_phonons, include_photons, which_ee_scattering, which_eQ_scattering, &
    &   max_lgrho_table, min_lgrho_table)
        use iso_fortran_env, only: error_unit
        implicit none
        integer, intent(in) :: handle
        logical, intent(in), optional :: include_electrons
        logical, intent(in), optional :: include_neutrons
        logical, intent(in), optional :: include_superfluid_phonons
        logical, intent(in), optional :: include_photons
        integer, intent(in), optional :: which_ee_scattering
        integer, intent(in), optional :: which_eQ_scattering
        real(dp), intent(in), optional :: max_lgrho_table
        real(dp), intent(in), optional :: min_lgrho_table
        type(conductivity_general_info), pointer :: rq
        integer :: ierr
        
        call get_conductivity_ptr(handle,rq,ierr)
        if (ierr /= 0) then
            write(error_unit,*) 'conductivity_set_controls: bad handle'
            return
        end if
        if (.not. rq% in_use) then
            write(error_unit,*) 'conductivity_set_controls: handle is not in use'
            return
        end if
        
        if (present(include_electrons)) rq% include_electrons = include_electrons
        if (present(include_neutrons)) rq% include_neutrons = include_neutrons
        if (present(include_superfluid_phonons)) rq% include_superfluid_phonons = &
        &   include_superfluid_phonons
        if (present(include_photons)) rq% include_photons = include_photons
        if (present(which_ee_scattering)) then
            if (which_ee_scattering==icond_sy06 .or. which_ee_scattering==icond_pcy) then
                rq% ee_scattering_fmla = which_ee_scattering
            else
                write(error_unit,'(a,i0,a)') 'unknown flag ',which_ee_scattering,' for ee scattering'
            end if
        end if
        if (present(which_eQ_scattering)) then
            if (which_eQ_scattering==icond_eQ_potekhin .or. which_eQ_scattering==icond_eQ_page) then
                rq% eQ_scattering_fmla = which_eQ_scattering
            else
                write(error_unit,'(a,i0,a)') 'unknown flag ',which_eQ_scattering,' for eQ scattering'
            end if
        end if
        if (present(max_lgrho_table)) rq% tab_off_lgrho = max_lgrho_table
        if (present(min_lgrho_table)) rq% tab_on_lgrho = min_lgrho_table
    end subroutine conductivity_set_controls

    subroutine get_thermal_conductivity( &
    &   rho,T,chi,Gamma,eta,mu_e,ionic,Tcn,K,which_components, &
    &   use_pcy,use_page)
        use nucchem_def, only: composition_info_type
        use eval_conductivity
        real(dp), intent(in) :: rho,T,Gamma,eta, mu_e, chi
        type(composition_info_type), intent(in) :: ionic
        real(dp), intent(in) :: Tcn ! neutron critical temperature
        type(conductivity_components), intent(out) :: K
        logical, intent(in), optional :: use_pcy
        logical, intent(in), optional :: use_page 
        logical, intent(in), optional :: &
        &   which_components(num_conductivity_channels)
        integer :: which_ee, which_eQ
        logical, dimension(num_conductivity_channels) :: K_components

        which_ee = icond_sy06
        which_eQ = icond_eQ_potekhin
        K_components = cond_use_all

        if (present(use_pcy)) then
            if (use_pcy) which_ee = icond_pcy
        end if

        if (present(use_page)) then
            if (use_page) which_eQ = icond_eQ_page
        end if

        if (present(which_components)) K_components = which_components

        call conductivity(rho,T,chi,Gamma,eta,mu_e,ionic,Tcn,K, &
        &   which_ee,which_eQ,K_components)
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
