    logical :: load_model_file
    character(len=64) :: model_file
    
    integer :: write_interval_for_terminal
    integer :: write_interval_for_terminal_header
    integer :: write_interval_for_history
    integer :: write_interval_for_profile
    integer :: starting_number_for_profile
    
    character(len=256) :: output_directory
   
    ! integration controls
    character(len=32) :: which_solver
    integer :: maximum_number_of_models
    real(dp) :: maximum_timestep
    real(dp) :: maximum_end_time
    real(dp) :: integration_tolerance
    real(dp) :: min_lg_temperature_integration
    real(dp) :: max_lg_temperature_integration
    
    real(dp), dimension(num_extra_real_controls) :: extra_real_controls
    integer, dimension(num_extra_integer_controls) :: extra_integer_controls
    logical, dimension(num_extra_logical_controls) :: extra_logical_controls
    
    logical :: fix_core_temperature
    real(dp) :: core_temperature
    
    logical :: fix_atmosphere_temperature_when_accreting
    real(dp) :: atmosphere_temperature_when_accreting
    real(dp) :: Mdot
    
    real(dp) :: start_time
    
    real(dp) :: core_mass
    real(dp) :: core_radius
    
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
    
    !  controls the shell Urca cooling
    logical :: turn_on_shell_Urca
    real(dp) :: shell_Urca_luminosity_coeff
    real(dp) :: lgP_shell_Urca
    
    ! controls for the atmosphere
    real(dp) :: lg_atm_light_element_column  ! lg(cgs), for the PCY atmosphere fit
    
    ! controls for the composition
    logical :: fix_Qimp
    real(dp) :: Qimp
