program test_NScool
    use iso_fortran_env, only : output_unit, error_unit
    use NScool_def
    use NScool_lib
    use argparse

    character(len=*), parameter :: default_dStar_dir = '../../'
    character(len=*), parameter :: default_inlist_file = 'test_inlist'
    character(len=64) :: my_dStar_dir
    character(len=64) :: inlist
    type(NScool_info), pointer :: s=>null()
    integer :: ierr, NScool_id, i
    
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
    ierr = 0

    call NScool_create_model(NScool_id,ierr)
    call check_okay('NScool_create_model',ierr)
    
    call NScool_evolve_model(NScool_id,ierr)
    call check_okay('NScool_evolve_model',ierr)

    ! write out the observed effective temperature at monitoring points
    ! (end of the epochs that were specified in the inlist).
    call get_NScool_info_ptr(NScool_id,s,ierr)    
    do i = 1, s% number_epochs
        write(output_unit,'(3(a," = ",es11.4))')  &
        &   't/d',s% t_monitor(i),'; obs. Teff/K', s% Teff_monitor(i), &
        & '; local Qb', s% Qb_monitor(i)
    end do
    
    call NScool_shutdown    
    write (output_unit,'(/,/,a,i3)') 'test_NScool exited with ierr = ',ierr
    
contains
	subroutine check_okay(msg,ierr)
		character(len=*), intent(in) :: msg
		integer, intent(inout) :: ierr
		if (ierr /= 0) then
			write (error_unit,*) trim(msg)//': ierr = ',ierr
			if (ierr < 0) stop
        else
            write (error_unit,*) trim(msg)//': okay'
		end if
	end subroutine check_okay
end program test_NScool
