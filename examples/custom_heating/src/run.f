program run_dStar
    use iso_fortran_env, only: output_unit
	use exceptions_lib
    use NScool_def
    use NScool_lib
    use alt_micro
    use argparse
    
    character(len=*), parameter :: default_inlist_file = 'inlist'
    character(len=64) :: my_dStar_dir
    character(len=64) :: inlist
    type(NScool_info), pointer :: s=>null()
    integer :: ierr, NScool_id, i
    real(dp), dimension(11) :: pred_Teff
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
    
    call get_NScool_info_ptr(NScool_id,s,ierr)
    if (check_okay% raised(ierr,'get_NScool_info_ptr')) stop
 
    s% other_set_Qimp => alt_Qimp
    s% other_set_heating => alt_heating

    call NScool_create_model(NScool_id,ierr)
    if (check_okay% raised(ierr,'NScool_create_model')) stop

    call NScool_evolve_model(NScool_id,ierr)        
    if (check_okay% raised(ierr,'NScool_evolve_model')) stop
   
    ! we don't want to compare the effective temp. at t = 0, the end of the 
    ! outburst, so we'll use indices 2:-
    pred_Teff = s% Teff_monitor(2:)/1.0e6    
    write(output_unit,*)
    write(output_unit,'(a7,a6)') 'time','Teff'
    write(output_unit,'(a7,a6)') '[d]','[MK]'
    write(output_unit,'(13("-"))')
    do i = 1, 11
        write(output_unit,'(f7.1,f6.3)')  &
        & s% t_monitor(i+1), pred_Teff(i)
    end do
    
    call NScool_shutdown
    
end program run_dStar
