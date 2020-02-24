program run_dStar
	use iso_fortran_env, only : error_unit
    use NScool_def
    use NScool_lib
    use argparse
    
    character(len=*), parameter :: default_dStar_dir = '../../dStar'
    character(len=*), parameter :: default_inlist_file = 'inlist'
    character(len=64) :: my_dStar_dir
    character(len=64) :: inlist
    type(NScool_info), pointer :: s=>null()
    integer :: ierr, NScool_id
    
    ierr = 0
    
    call command_arg_set( &
        & 'dStar_directory',"sets the main dStar root directory",ierr, &
        & flag='D',takes_parameter=.TRUE.)
    call check_okay('set command argument dStar_directory',ierr)
    
    call command_arg_set( &
        & 'inlist_file','sets the namelist parameter file',ierr, &
        & flag='I',takes_parameter=.TRUE.)
    call check_okay('set command argument inlist file',ierr)
    
    call parse_arguments(ierr)
    call check_okay('parse_arguments',ierr)

    my_dStar_dir = trim(get_command_arg('dStar_directory'))
    if (len_trim(my_dStar_dir)==0) my_dStar_dir = default_dStar_dir
    inlist = trim(get_command_arg('inlist_file'))
    if (len_trim(inlist)==0) inlist = default_inlist_file
    
    call NScool_init(my_dStar_dir, ierr)
    call check_okay('NScool_init',ierr)
    
    NScool_id = alloc_NScool(ierr)
    call check_okay('NScool_id',ierr)
    
    call NScool_setup(NScool_id,inlist,ierr)
    call check_okay('NScool_setup',ierr)
    
    call NScool_create_model(NScool_id,ierr)
    call check_okay('NScool_create_model',ierr)

    call NScool_evolve_model(NScool_id,ierr)        
    call check_okay('NScool_evolve_model',ierr)

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
