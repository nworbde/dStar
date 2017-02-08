program run_dStar
	use iso_fortran_env, only : output_unit, error_unit
    use NScool_def
    use NScool_lib
    use constants_def, only : boltzmann
    use argparse
    
    character(len=*), parameter :: default_dStar_dir = '../../dStar'
    character(len=*), parameter :: default_inlist_file = 'inlist'
    character(len=*), parameter :: default_Q_heating_shallow = '0.0' ! MeV
    character(len=*), parameter :: default_core_mass = '1.6' ! Msun
    character(len=*), parameter :: default_core_radius = '11.0' ! km
    
    character(len=64) :: my_dStar_dir, inlist, Q_heating_shallow_arg, &
        & core_mass_arg, core_radius_arg
    
    real(dp) :: Q_heating_shallow, core_mass, core_radius
    type(NScool_info), pointer :: s=>null()
    integer :: ierr, NScool_id, i
    
    ierr = 0
    call command_arg_set( &
        & 'dStar_directory',"sets the main dStar root directory",ierr, &
        & flag='D',takes_parameter=.TRUE.)
    call check_okay('set command argument dStar_directory',ierr)
    
    call command_arg_set( &
        & 'inlist_file','sets the namelist parameter file',ierr, &
        & flag='I',takes_parameter=.TRUE.)
    call check_okay('set command argument inlist file',ierr)
    
    call command_arg_set( &
        & 'Qheating','shallow heating per nucleon',ierr, &
        & flag='Q',takes_parameter=.TRUE.)
    call check_okay('set command argument shallow heating',ierr)

    call command_arg_set( &
        & 'Mcore','core mass',ierr, &
        & flag='M',takes_parameter=.TRUE.)
    call check_okay('set command argument core mass',ierr)

    call command_arg_set( &
        & 'Rcore','core radius',ierr, &
        & flag='R',takes_parameter=.TRUE.)
    call check_okay('set command argument core radius',ierr)
    
    call parse_arguments(ierr)
    call check_okay('parse_arguments',ierr)

    my_dStar_dir = trim(get_command_arg('dStar_directory'))
    if (len_trim(my_dStar_dir)==0) my_dStar_dir = default_dStar_dir

    inlist = trim(get_command_arg('inlist_file'))
    if (len_trim(inlist)==0) inlist = default_inlist_file

    Q_heating_shallow_arg = trim(get_command_arg('Qheating'))
    if (len_trim(Q_heating_shallow_arg)==0)  &
        & Q_heating_shallow_arg = default_Q_heating_shallow
    read(Q_heating_shallow_arg,*) Q_heating_shallow

    core_mass_arg = trim(get_command_arg('Mcore'))
    if (len_trim(core_mass_arg)==0)  &
        & core_mass_arg = default_core_mass
    read(core_mass_arg,*) core_mass
    print *, 'core mass = ',core_mass

    core_radius_arg = trim(get_command_arg('Rcore'))
    if (len_trim(core_radius_arg)==0)  &
        & core_radius_arg = default_core_radius
    read(core_radius_arg,*) core_radius
    
    call NScool_init(my_dStar_dir, ierr)
    call check_okay('NScool_init',ierr)
    
    NScool_id = alloc_NScool(ierr)
    call check_okay('NScool_id',ierr)

    call NScool_setup(NScool_id,inlist,ierr)
    call check_okay('NScool_setup',ierr)
    
    call get_NScool_info_ptr(NScool_id,s,ierr)
    call check_okay('get_NScool_info_ptr',ierr)
    s% Q_heating_shallow = Q_heating_shallow
    s% core_mass = core_mass
    s$ core_radius = core_radius
    s% Mcore = core_mass
    s% Rcore = core_radius
    
    call NScool_create_model(NScool_id,ierr)
    call check_okay('NScool_create_model',ierr)
    
    call NScool_evolve_model(NScool_id,ierr)        
    call check_okay('NScool_evolve_model',ierr)

    write(output_unit,*)
    write(output_unit,'(a7,a6,a11)') 'time','Teff','Flux/mdot'
    write(output_unit,'(a7,a6,a11)') '[d]','[MK]','[MeV/u]'
    write(output_unit,'(24("-"))')
    do i = 1, s% number_epochs
        write(output_unit,'(f7.1,f6.3,es11.4)')  &
        & s% t_monitor(i), s% Teff_monitor(i)/1.0e6, s% Qb_monitor(i)
    end do
    
    call NScool_shutdown
    
contains
	subroutine check_okay(msg,ierr)
		character(len=*), intent(in) :: msg
		integer, intent(inout) :: ierr
		if (ierr /= 0) then
			write (error_unit,*) trim(msg)//': ierr = ',ierr
			if (ierr < 0) stop
		end if
	end subroutine check_okay
end program run_dStar
