program run_dStar
	use exceptions_lib
    use NScool_def
    use NScool_lib
    use argparse
    
    character(len=*), parameter :: default_inlist_file = 'inlist'
    character(len=64) :: my_dStar_dir
    character(len=64) :: inlist
    type(NScool_info), pointer :: s=>null()
    integer :: ierr, NScool_id
    type(failure) :: check_okay=failure(scope='run_dStar')
    
    ierr = 0
    
    call command_arg_set( &
        & 'dStar_directory',"sets the main dStar root directory",ierr, &
        & flag='D',takes_parameter=.TRUE.)
    if (check_okay% raised(ierr,'set command argument dStar_directory')) stop
    
    call command_arg_set( &
        & 'inlist_file','sets the namelist parameter file',ierr, &
        & flag='I',takes_parameter=.TRUE.)
    if (check_okay% raised(ierr,'set command argument inlist file')) stop
    
    call parse_arguments(ierr)
    if (check_okay% raised(ierr,'parse_arguments')) stop

    my_dStar_dir = trim(get_command_arg('dStar_directory'))
    inlist = trim(get_command_arg('inlist_file'))
    if (len_trim(inlist)==0) inlist = default_inlist_file
    
    call NScool_init(my_dStar_dir, ierr)
    if (check_okay% raised(ierr,'NScool_init')) stop
    
    NScool_id = alloc_NScool(ierr)
    if (check_okay% raised(ierr,'NScool_id')) stop
    
    call NScool_setup(NScool_id,inlist,ierr)
    if (check_okay% raised(ierr,'NScool_setup')) stop
    
    call NScool_create_model(NScool_id,ierr)
    if (check_okay% raised(ierr,'NScool_create_model')) stop

    call NScool_evolve_model(NScool_id,ierr)        
    if (check_okay% raised(ierr,'NScool_evolve_model')) stop

    call NScool_shutdown
    
end program run_dStar
