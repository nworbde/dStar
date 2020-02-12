module NScool_ctrls_io
    use NScool_def
    include 'NScool_controls.inc'

    namelist /controls/ &
        load_model_file, & 
        model_file, & 
        write_interval_for_terminal, & 
        write_interval_for_terminal_header, & 
        write_interval_for_history, & 
        write_interval_for_profile, & 
        starting_number_for_profile, &
        suppress_first_step_output, &
        output_directory, & 
        which_solver, & 
        maximum_number_of_models, & 
        maximum_timestep, & 
        integration_tolerance, & 
        min_lg_temperature_integration, &
        max_lg_temperature_integration, &
        use_other_set_Qimp, &
        use_other_sf_critical_temperatures, &
        extra_real_controls, & 
        extra_integer_controls, & 
        extra_logical_controls, & 
        fix_core_temperature, & 
        core_temperature, & 
        make_inner_boundary_insulating, &
        fix_atmosphere_temperature_when_accreting, & 
        atmosphere_temperature_when_accreting, &
        number_epochs, &
        epoch_Mdots, &
        epoch_boundaries, &
        core_mass, &
        core_radius, &
        lgPcrust_bot, &
        lgPcrust_top, &
        target_resolution_lnP,  &
        lgP_min_heating_outer, &
        lgP_max_heating_outer, &
        Q_heating_outer, &
        lgP_min_heating_inner, &
        lgP_max_heating_inner, &
        Q_heating_inner, &
        turn_on_extra_heating, &
        lgP_min_heating_shallow, &
        lgP_max_heating_shallow, &
        Q_heating_shallow, &
        use_other_set_heating, &
        which_proton_1S0_gap, &
        which_neutron_1S0_gap, &
        which_neutron_3P2_gap, &
        scale_sf_critical_temperatures, &
        use_pcy_for_ee_scattering, &
        use_page_for_eQ_scattering, &
        use_ee_conductivity, &
        use_ei_conductivity, &
        use_eQ_conductivity, &
        use_sf_conductivity, &
        use_rad_opacity, &
        eos_gamma_melt_pt, &
        eos_rsi_melt_pt, &
        eos_nuclide_abundance_threshold, &
        eos_pasta_transition_in_fm3, &
        eos_cluster_transition_in_fm3, &
        use_core_nu_bremsstrahlung, &
        use_core_nu_mUrca, &
        use_core_nu_dUrca, &
        use_core_nu_PBF, &
        use_crust_nu_pair, &
        use_crust_nu_photo, &
        use_crust_nu_plasma, &
        use_crust_nu_bremsstrahlung, &
        use_crust_nu_pbf, &
        turn_on_shell_Urca, &
        shell_Urca_luminosity_coeff, &
        lgP_shell_Urca, &
        lg_atm_light_element_column, &
        atm_model, &
        fix_Qimp, &
        Qimp
        
contains
    subroutine do_one_setup(id,inlist,ierr)
        integer, intent(in) :: id
        character(len=*), intent(in) :: inlist
        integer, intent(out) :: ierr
        type(NScool_info), pointer :: s
        
        call get_NScool_info_ptr(id,s,ierr)
        if (ierr /= 0) return
        call set_default_controls
        call read_controls(id,inlist,ierr)
        call store_controls(s, ierr)
    end subroutine do_one_setup

    subroutine read_controls(id,filename,ierr)
        integer, intent(in) :: id
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ierr
        type(NScool_info), pointer :: s
        integer :: iounit

        ierr = 0
        call get_NScool_info_ptr(id,s,ierr)
        if (ierr /= 0) return

        if (len_trim(filename) > 0) then
            open (newunit=iounit,file=trim(filename),  &
            & action='read', delim='quote', status='old', iostat=ierr)
            if (ierr /= 0) then
                write(*,*) 'failed to open control namelist file ',trim(filename)
                return
            end if
            read(iounit,nml=controls,iostat=ierr)
            close(iounit)
            if (ierr /= 0) then
                write(*,'(///,a)') 'failed while reading control namelist file '//trim(filename)

                ! the following will generate a more informative description of the failure
                close(iounit)
                open (newunit=iounit,file=trim(filename),  &
                & action='read', delim='quote', status='old', iostat=ierr)
                read(iounit,nml=controls)
                close(iounit)
                stop
            end if
        end if
    end subroutine read_controls
   
    subroutine set_default_controls()
        include 'controls.defaults'
    end subroutine set_default_controls

    subroutine store_controls(s,ierr)
        use constants_def, only : Msun, julian_day
        use num_lib, only : solver_option
        type(NScool_info), pointer :: s
        integer, intent(out) :: ierr

        ierr = 0

        s% load_model_file = load_model_file
        s% model_file = model_file
        s% write_interval_for_terminal = write_interval_for_terminal
        s% write_interval_for_terminal_header = write_interval_for_terminal_header
        s% write_interval_for_history = write_interval_for_history
        s% write_interval_for_profile = write_interval_for_profile
        s% starting_number_for_profile = starting_number_for_profile
        s% suppress_first_step_output = suppress_first_step_output
        s% output_directory = output_directory
        s% which_solver = trim(which_solver)
        s% base_profile_filename = trim(s% output_directory)//'/profile'
        s% history_filename = trim(s% output_directory)//'/history.data'
        s% profile_manifest_filename = trim(s% output_directory)//'/profiles'

        s% maximum_number_of_models = maximum_number_of_models
        s% maximum_timestep = maximum_timestep*julian_day
        s% integration_tolerance = integration_tolerance
        s% min_lg_temperature_integration = min_lg_temperature_integration
        s% max_lg_temperature_integration = max_lg_temperature_integration

        s% use_other_set_Qimp = use_other_set_Qimp
        s% use_other_sf_critical_temperatures = use_other_sf_critical_temperatures
        s% extra_real_controls = extra_real_controls
        s% extra_integer_controls = extra_integer_controls
        s% extra_logical_controls = extra_logical_controls

        s% fix_core_temperature = fix_core_temperature
        s% core_temperature = core_temperature
        s% make_inner_boundary_insulating = make_inner_boundary_insulating
        s% fix_atmosphere_temperature_when_accreting = fix_atmosphere_temperature_when_accreting
        s% atmosphere_temperature_when_accreting = atmosphere_temperature_when_accreting
        s% epoch_data = epoch_data
!         s% number_epochs = number_epochs
!         s% epoch_Mdots = epoch_Mdots
!         s% epoch_boundaries = epoch_boundaries*julian_day

        s% Mcore = core_mass
        s% Rcore = core_radius
        s% Psurf = 10.0**lgPcrust_top
        s% Pcore = 10.0**lgPcrust_bot
        s% target_resolution_lnP = target_resolution_lnP
        s% Tcore = core_temperature
        s% lgP_min_heating_outer = lgP_min_heating_outer
        s% lgP_max_heating_outer = lgP_max_heating_outer
        s% Q_heating_outer = Q_heating_outer
        s% lgP_min_heating_inner = lgP_min_heating_inner
        s% lgP_max_heating_inner = lgP_max_heating_inner
        s% Q_heating_inner = Q_heating_inner
        s% turn_on_extra_heating = turn_on_extra_heating
        s% lgP_min_heating_shallow = lgP_min_heating_shallow
        s% lgP_max_heating_shallow = lgP_max_heating_shallow
        s% Q_heating_shallow = Q_heating_shallow
        s% use_other_set_heating = use_other_set_heating

        s% which_proton_1S0_gap = which_proton_1S0_gap
        s% which_neutron_1S0_gap = which_neutron_1S0_gap
        s% which_neutron_3P2_gap = which_neutron_3P2_gap
        s% scale_sf_critical_temperatures = scale_sf_critical_temperatures

        s% use_pcy_for_ee_scattering = use_pcy_for_ee_scattering
        s% use_page_for_eQ_scattering = use_page_for_eQ_scattering

        s% use_ee_conductivity = use_ee_conductivity
        s% use_ei_conductivity = use_ei_conductivity
        s% use_eQ_conductivity = use_eQ_conductivity
        s% use_sf_conductivity = use_sf_conductivity
        s% use_rad_opacity = use_rad_opacity

        s% eos_gamma_melt_pt = eos_gamma_melt_pt
        s% eos_rsi_melt_pt = eos_rsi_melt_pt
        s% eos_nuclide_abundance_threshold = eos_nuclide_abundance_threshold
        s% eos_pasta_transition_in_fm3 = eos_pasta_transition_in_fm3
        s% eos_cluster_transition_in_fm3 = eos_cluster_transition_in_fm3

        s% use_core_nu_bremsstrahlung = use_core_nu_bremsstrahlung
        s% use_core_nu_mUrca = use_core_nu_mUrca
        s% use_core_nu_dUrca = use_core_nu_dUrca
        s% use_core_nu_PBF = use_core_nu_PBF
        s% use_crust_nu_pair = use_crust_nu_pair
        s% use_crust_nu_photo = use_crust_nu_photo
        s% use_crust_nu_plasma = use_crust_nu_plasma
        s% use_crust_nu_bremsstrahlung = use_crust_nu_bremsstrahlung
        s% use_crust_nu_pbf = use_crust_nu_pbf

        s% turn_on_shell_Urca = turn_on_shell_Urca
        s% shell_Urca_luminosity_coeff = shell_Urca_luminosity_coeff
        s% lgP_shell_Urca = lgP_shell_Urca

        s% lg_atm_light_element_column = lg_atm_light_element_column
        s% atm_model = atm_model
        s% fix_Qimp = fix_Qimp
        s% Qimp = Qimp

   end subroutine store_controls
   
   subroutine write_controls(io,ierr)
       integer, intent(in) :: io
       integer, intent(out) :: ierr
       write(io,nml=controls,iostat=ierr)
   end subroutine write_controls
   
end module NScool_ctrls_io
