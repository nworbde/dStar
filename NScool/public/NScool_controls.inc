    logical :: load_model_file
    character(len=64) :: model_file
    
    integer :: write_interval_for_terminal
    integer :: write_interval_for_terminal_header
    integer :: write_interval_for_history
    integer :: write_interval_for_profile
    integer :: starting_number_for_profile
    logical :: suppress_first_step_output
    
    character(len=256) :: output_directory
   
    ! integration controls
    character(len=32) :: which_solver
    integer :: maximum_number_of_models
    real(dp) :: maximum_timestep
    real(dp) :: integration_tolerance
    real(dp) :: min_lg_temperature_integration
    real(dp) :: max_lg_temperature_integration
    
    logical :: fix_core_temperature
    real(dp) :: core_temperature
    logical :: make_inner_boundary_insulating
    
    logical :: fix_atmosphere_temperature_when_accreting
    real(dp) :: atmosphere_temperature_when_accreting
    
    integer :: number_epochs
    real(dp), dimension(max_number_epochs) :: epoch_Mdots
    real(dp), dimension(0:max_number_epochs) :: epoch_boundaries
    
    real(dp) :: core_mass
    real(dp) :: core_radius
    
    ! lgPcrust_top is also used for the bottom of the atmosphere with the bc09 model
    real(dp) :: lgPcrust_bot
    real(dp) :: lgPcrust_top
    
    real(dp) :: target_resolution_lnP
    
    ! controls for the crust heating
    real(dp) :: lgP_min_heating_outer
    real(dp) :: lgP_max_heating_outer
    real(dp) :: Q_heating_outer ! MeV
    real(dp) :: lgP_min_heating_inner
    real(dp) :: lgP_max_heating_inner
    real(dp) :: Q_heating_inner ! MeV
    logical  :: turn_on_extra_heating
    real(dp) :: lgP_min_heating_shallow
    real(dp) :: lgP_max_heating_shallow
    real(dp) :: Q_heating_shallow ! MeV
    
    ! controls for the superfluidity
    logical :: use_other_sf_critical_temperatures
    character(len=32) :: which_proton_1S0_gap
    character(len=32) :: which_neutron_1S0_gap
    character(len=32) :: which_neutron_3P2_gap
    real(dp), dimension(3) :: scale_sf_critical_temperatures
    
    ! controls for the conductivity: setting these = .TRUE. means don't use the defaults
    logical :: use_pcy_for_ee_scattering
    logical :: use_page_for_eQ_scattering
    ! set which channels to include
    logical :: use_ee_conductivity
    logical :: use_ei_conductivity
    logical :: use_eQ_conductivity
    logical :: use_sf_conductivity
    logical :: use_rad_opacity
    
    ! EOS controls: set to < 0 to use defaults
    real(dp) :: eos_gamma_melt_pt
    real(dp) :: eos_rsi_melt_pt
    real(dp) :: eos_nuclide_abundance_threshold
    real(dp) :: eos_pasta_transition_in_fm3
    real(dp) :: eos_cluster_transition_in_fm3
    
    ! neutrino controls: set which channels to include
    logical :: use_core_nu_bremsstrahlung
    logical :: use_core_nu_mUrca
    logical :: use_core_nu_dUrca
    logical :: use_core_nu_PBF
    logical :: use_crust_nu_pair
    logical :: use_crust_nu_photo
    logical :: use_crust_nu_plasma
    logical :: use_crust_nu_bremsstrahlung
    logical :: use_crust_nu_pbf
    
    !  controls the shell Urca cooling
    logical :: turn_on_shell_Urca
    real(dp) :: shell_Urca_luminosity_coeff
    real(dp) :: lgP_shell_Urca
    
    ! controls for the atmosphere
    real(dp) :: lg_atm_light_element_column  ! lg(cgs), for the atmosphere fit
    character(len=8) :: atm_model
    
    ! controls for the composition
    logical :: use_other_set_Qimp
    logical :: fix_Qimp
    real(dp) :: Qimp

    ! controls for user-defined routines
    real(dp), dimension(num_extra_real_controls) :: extra_real_controls
    integer, dimension(num_extra_integer_controls) :: extra_integer_controls
    logical, dimension(num_extra_logical_controls) :: extra_logical_controls
