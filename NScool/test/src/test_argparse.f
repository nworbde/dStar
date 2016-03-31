program test_argparse
    use argparse
    integer :: i, ierr
    call command_arg_set( &
        & 'dStar_directory',"sets the main dStar root directory",ierr, &
        & flag='D',takes_parameter=.TRUE.)
    if (failure('command_arg_set dStar_directory')) stop

    call command_arg_set( &
        & 'setQimp',"sets the Qimp flag",ierr,flag='Q',takes_parameter=.FALSE.)
    if (failure('command_arg_set setQimp')) stop

    call command_arg_set( &
        & 'inlist',"sets the inlist",ierr)
    if (failure('command_arg_set inlist')) stop

    call command_arg_set( &
        & 'other_inlist',"sets the other inlist",ierr)
    if (failure('command_arg_set inlist')) stop
    
    call parse_arguments(ierr)
    if (failure('parse_arguments')) stop
    
    write(*,*) 'dStar_directory = ',get_command_arg('dStar_directory')
    write(*,*) 'setQimp = ',get_command_arg('setQimp')
    write(*,*) 'inlist =',get_command_arg('inlist')
    write(*,*) 'other inlist =',get_command_arg('other_inlist')
    
contains
    function failure(msg)
        character(len=*), intent(in) :: msg
        logical :: failure
        
        failure = .FALSE.
        if (ierr /= 0) then
            failure = .TRUE.
            write(*,*) 'failure in '//trim(msg)
        end if
    end function failure
end program test_argparse
