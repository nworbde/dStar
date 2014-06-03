    logical :: load_model_file
    character(len=64) :: model file
    
    integer :: write_interval_for_terminal
    integer :: write_interval_for_terminal_header
    integer :: write_interval_for_history
    integer :: write_interval_for_profile
    
    character(len=256) :: output_directory
   
    character(len=32) :: which_solver
    integer :: maximum_number_of_models
    real(dp) :: maximum_timestep
    real(dp) :: maximum_end_time
    