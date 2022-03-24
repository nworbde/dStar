program run_dStar
	use exceptions_lib
    use NScool_def
    use NScool_lib
    use constants_def, only : boltzmann
    use argparse
    
    character(len=*), parameter :: default_dStar_dir = '../../dStar'
    character(len=*), parameter :: default_inlist_file = 'inlist'
    real(dp), parameter :: default_Q_heating_shallow = 0.0 ! MeV
    real(dp), parameter :: default_core_mass = 1.6 ! Msun
    real(dp), parameter :: default_core_radius = 11.0 ! km
    
    character(len=64) :: my_dStar_dir, inlist, Q_heating_shallow_arg, &
        & core_mass_arg, core_radius_arg
    
    real(dp) :: Q_heating_shallow, core_mass, core_radius
    type(NScool_info), pointer :: s=>null()
    integer :: ierr, NScool_id, i
    type(failure) :: check_okay=failure(scope='run_dStar')
    
    ierr = 0
    call command_arg_set( &
        & 'dStar_directory',"sets the main dStar root directory",ierr, &
        & flag='D',takes_parameter=.TRUE.)
    if (check_okay% failure(ierr,'set command argument dStar_directory')) stop
    
    call command_arg_set( &
        & 'inlist_file','sets the namelist parameter file',ierr, &
        & flag='I',takes_parameter=.TRUE.)
    if (check_okay% failure(ierr,'set command argument inlist file')) stop
    
    call command_arg_set( &
        & 'Qheating','shallow heating per nucleon',ierr, &
        & flag='Q',takes_parameter=.TRUE.)
    if (check_okay% failure(ierr,'set command argument shallow heating')) stop

    call command_arg_set( &
        & 'Mcore','core mass',ierr, &
        & flag='M',takes_parameter=.TRUE.)
    if (check_okay% failure(ierr,'set command argument core mass')) stop

    call command_arg_set( &
        & 'Rcore','core radius',ierr, &
        & flag='R',takes_parameter=.TRUE.)
    if (check_okay% failure(ierr,'set command argument core radius')) stop
    
    call parse_arguments(ierr)
    if (check_okay% failure(ierr,'parse_arguments')) stop

    my_dStar_dir = trim(get_command_arg('dStar_directory'))
    inlist = trim(get_command_arg('inlist_file'))
    if (len_trim(inlist)==0) inlist = default_inlist_file

    Q_heating_shallow_arg = trim(get_command_arg('Qheating'))
    if (len_trim(Q_heating_shallow_arg) > 0)  then
        read(Q_heating_shallow_arg,*) Q_heating_shallow
    else
        Q_heating_shallow = default_Q_heating_shallow
    end if

    core_mass_arg = trim(get_command_arg('Mcore'))
    if (len_trim(core_mass_arg) > 0)  then
        read(core_mass_arg,*) core_mass
    else
        core_mass = default_core_mass
    end if
    
    core_radius_arg = trim(get_command_arg('Rcore'))
    if (len_trim(core_radius_arg) > 0)  then
        read(core_radius_arg,*) core_radius
    else
        core_radius = default_core_radius
    end if
    
    call NScool_init(my_dStar_dir, ierr)
    if (check_okay% failure(ierr,'NScool_init')) stop
    
    NScool_id = alloc_NScool(ierr)
    if (check_okay% failure(ierr,'NScool_id')) stop

    call NScool_setup(NScool_id,inlist,ierr)
    if (check_okay% failure(ierr,'NScool_setup')) stop
    
    call get_NScool_info_ptr(NScool_id,s,ierr)
    if (check_okay% failure(ierr,'get_NScool_info_ptr')) stop

    s% Q_heating_shallow = Q_heating_shallow
    s% core_mass = core_mass
    s% core_radius = core_radius
    s% Mcore = core_mass
    s% Rcore = core_radius
    
    call NScool_create_model(NScool_id,ierr)
    if (check_okay% failure(ierr,'NScool_create_model')) stop
    
    call NScool_evolve_model(NScool_id,ierr)        
    if (check_okay% failure(ierr,'NScool_evolve_model')) stop

    write(output_unit,*)
    write(output_unit,'(a7,a6,a11)') 'time','Teff','Flux/mdot'
    write(output_unit,'(a7,a6,a11)') '[d]','[MK]','[MeV/u]'
    write(output_unit,'(24("-"))')
    do i = 1, s% number_epochs
        write(output_unit,'(f7.1,f6.3,es11.4)')  &
        & s% t_monitor(i), s% Teff_monitor(i)/1.0e6, s% Qb_monitor(i)
    end do
    
    call NScool_shutdown
end program run_dStar
