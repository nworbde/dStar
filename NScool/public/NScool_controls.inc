    logical :: load_model_file
    character(len=64) :: model_file
    
    integer :: write_interval_for_terminal
    integer :: write_interval_for_terminal_header
    integer :: write_interval_for_history
    integer :: write_interval_for_profile
    
    character(len=256) :: output_directory
   
    character(len=32) :: which_solver
    integer :: maximum_number_of_models
    real(dp) :: maximum_timestep
    real(dp) :: maximum_end_time
    
    real(dp), dimension(num_extra_real_controls) :: extra_real_controls
    integer, dimension(num_extra_integer_controls) :: extra_integer_controls
    logical, dimension(num_extra_logical_controls) :: extra_logical_controls
    
    logical :: fix_core_temperature
    real(dp) :: core_temperature
    
    logical :: fix_atmosphere_temperature_when_accreting
    real(dp) :: atmosphere_temperature_when_accreting
    