module argparse
    use utils_def
    use utils_lib
    
    integer, parameter :: max_command_args = 32
    integer, parameter :: doc_len = 256
    integer, parameter :: arg_len = 128
    
    type command_arg_type
        integer :: arg_id
        character(len=arg_len) :: argname
        character(len=doc_len) :: docstring
        character(len=arg_len) :: value
        integer :: position
        character(len=1) :: flag
        logical :: takes_parameter
        logical :: in_use
    end type command_arg_type

    character(len=arg_len), dimension(-3:-1), parameter :: parse_arguments_error_msgs =  &
    &   [ character(len=arg_len) :: 'wrong number positional parameters','unrecognized option', &
    &   'missing parameter']
    
    integer, save :: number_command_args = 0
    integer, save :: number_positional_args = 0
    integer, dimension(max_command_args), save :: positional_args
    type(command_arg_type), dimension(max_command_args), save :: command_args
    type(integer_dict), pointer, save :: arguments_dict=>null()
    type(integer_dict), pointer, save :: options_dict=>null()
    
contains
    
    subroutine command_arg_set( &
        & argname,doc,ierr,flag,takes_parameter)
        
        character(len=*), intent(in) :: argname
        character(len=*), intent(in) :: doc
        integer, intent(out) :: ierr
        character(len=1), optional :: flag
        logical, optional :: takes_parameter
        integer :: i
        integer :: current_command_arg
        logical :: positional

        ierr = 0
        current_command_arg = -1
        do i=1, max_command_args
            if (.not. command_args(i)% in_use) then
                command_args(i)% in_use = .TRUE.
                command_args(i)% arg_id = i
                current_command_arg = i
                exit
            end if
        end do
        if (current_command_arg == -1) then
            ierr = -1
            return
        end if
        if (command_args(current_command_arg)% arg_id /= current_command_arg) then
            ierr = -1
            return
        end if
        
        command_args(current_command_arg)% argname = argname
        command_args(current_command_arg)% docstring = doc
        command_args(current_command_arg)% value = ''
        command_args(current_command_arg)% flag = ''
        command_args(current_command_arg)% takes_parameter = .FALSE.

        if (present(flag))  &
        & command_args(current_command_arg)% flag = flag
        
        positional = (len_trim(command_args(current_command_arg)% flag) == 0)

        if (positional) then
            number_positional_args = number_positional_args + 1
            if (number_positional_args > max_command_args) then
                ierr = -2
                return
            end if
            command_args(current_command_arg)% position = number_positional_args
            positional_args(number_positional_args) = current_command_arg
        else
            command_args(current_command_arg)% flag = flag
            call integer_dict_define(options_dict,flag,current_command_arg,ierr)
            if (ierr /= 0) return
            command_args(current_command_arg)% position = -1
            if (present(takes_parameter)) &
            & command_args(current_command_arg)% takes_parameter = takes_parameter
            if (.not. command_args(current_command_arg)% takes_parameter)  &
            &  command_args(current_command_arg)% value = 'F'
        end if
        
        call integer_dict_define(arguments_dict, &
            & trim(adjustl(argname)),current_command_arg,ierr)
        
        number_command_args = current_command_arg
    end subroutine command_arg_set
    
    function get_command_arg(argname)
        character(len=*), intent(in) :: argname
        character(len=arg_len) :: get_command_arg
        integer :: indx, ierr
        
        call integer_dict_lookup(arguments_dict, argname, indx, ierr)
        if (ierr /= 0) then
            get_command_arg = ''
            return
        end if
        
        get_command_arg = command_args(indx)% value
    end function get_command_arg
    
    subroutine parse_arguments(ierr)
        use, intrinsic :: iso_fortran_env, only: error_unit
        integer, intent(out) :: ierr
        logical :: done_with_opt, show_help
        character(len=arg_len) :: tok
        character(len=1) :: op
        integer :: iop, ipos, num_args, indx
        ierr = 0

        iop = 1
        ipos = 1
        done_with_opt = .FALSE.
        show_help = .FALSE.
        
        do
            call get_command_argument(iop,tok)
            if (len_trim(tok) == 0) exit
            iop = iop + 1
            if (done_with_opt .or. .not. scan(tok,'-') == 1) then
                call store_pos_arg(ierr)
                if (ierr /= 0) then
                    ierr = -3
                    exit
                end if
            else
                op = tok(2:2)
                ! is it a known option?
                if (op == '-' .and. len_trim(tok)==2) then
                    done_with_opt = .TRUE.
                    cycle
                else if (op == 'h') then
                    show_help = .TRUE.
                    exit
                end if
                call integer_dict_lookup(options_dict, op, indx, ierr)
                if (ierr /= 0) then
                    ierr = -2
                    exit
                end if
                ! does it take a parameter?
                if (command_args(indx)% takes_parameter) then
                    if (len_trim(tok(3:)) > 0) then
                        command_args(indx)% value = tok(3:)
                    else
                        call get_command_argument(iop,tok)
                        if (len_trim(tok) == 0 .or. scan(tok,'-') == 1) then
                            ierr = -1
                            exit
                        end if
                        iop = iop+1
                        command_args(indx)% value = tok
                    end if
                else
                    command_args(indx)% value = 'T'
                end if
            end if

        end do

        ! test for error condition
        if (ierr /= 0) then
            write(error_unit,*) trim(parse_arguments_error_msgs(ierr))
        end if
        if (ierr /=0 .or. show_help) then
            call write_usage()
        end if

    contains
        subroutine store_pos_arg(ierr)
            integer, intent(out) :: ierr
            if (ipos <= number_command_args) then
                command_args(positional_args(ipos))% value = trim(tok)
                ipos = ipos +1
                ierr = 0
            else
                ierr  = -1
            end if
        end subroutine store_pos_arg
    end subroutine parse_arguments

    subroutine write_usage()
        use, intrinsic :: iso_fortran_env, only: error_unit
        character(len=doc_len) :: command_name
        character(len=max_command_args*number_command_args) :: command_line
        integer :: i
        
        call get_command_argument(0,command_name)
        command_line = "Usage: "//trim(command_name)
        do i = 1, number_command_args
            if (command_args(i)% flag /= '') then
                call append_op('-'//command_args(i)% flag)
                if (command_args(i)% takes_parameter) then
                    call append_op('<'//trim(command_args(i)% argname)//'>')
                end if
            else
                call append_op('<'//trim(command_args(i)% argname)//'>')
            end if
        end do
        write (error_unit,*)
        write(error_unit,*) trim(command_line)
        do i = 1, number_command_args
            if (command_args(i)% flag /= '') then
                if (command_args(i)% takes_parameter) then
                    write (error_unit,'(t4,"-",a," <",a,">:")') command_args(i)% flag,trim(command_args(i)% argname)
                else
                    write (error_unit,'(t4,"-",a,":")') command_args(i)% flag
                end if
            else
                write (error_unit,'(t4,"<",a,">:")') trim(command_args(i)% argname)
            end if
            write (error_unit,'(t8,a)') trim(command_args(i)% docstring)
        end do
        
    contains
        subroutine append_op(tok)
            character(len=*), intent(in) :: tok
            command_line = trim(command_line) // ' '//trim(tok)
        end subroutine append_op
    end subroutine write_usage

end module argparse
