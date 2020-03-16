module exceptions_lib
    use iso_fortran_env, only: error_unit
    implicit none
    
    integer, private, save :: verbosity=3
    integer, private, save :: warning_count = 0
    integer, private, save :: max_warnings = 24
    logical, private, save :: suppress_warnings_after_max = .FALSE.

    character(len=*), parameter, private :: indent4='(tr4,a)'
    character(len=*), parameter, private :: alert_fmt='(a,": ",a)'
    integer, parameter :: message_length = 256
    integer, parameter :: scope_length = 64
    
    type exception
        character(len=scope_length) :: scope = ''
        character(len=message_length) :: message = ''
    contains
        procedure :: set_scope
        procedure :: set_message
    end type exception

    type, extends(exception) :: failure
    contains
        procedure :: raised => failure_raise
    end type failure

    type, extends(exception) :: warning
    contains
        procedure :: raised => warning_raise
    end type warning
    
    type, extends(exception) :: alert
        integer :: level
    contains
        procedure :: set_level
        procedure :: report => alert_message
    end type alert
    
    type, extends(exception) :: assertion
        logical :: halt = .TRUE.
    contains
        procedure :: set_halt_on_assertion_failure
        procedure :: assert
    end type assertion
    
contains
    
    subroutine set_verbosity(v)
        integer, intent(in) :: v
        verbosity = v
    end subroutine set_verbosity
    
    subroutine set_max_warnings(max_count)
        integer, intent(in) :: max_count
        max_warnings = max_count
    end subroutine
    
    subroutine reset_warnings_count(count)
        integer, intent(in), optional :: count
        warning_count = 0
        if (present(count)) warning_count = count
    end subroutine reset_warnings_count
    
    subroutine set_suppress_warnings_after_max(suppress)
        logical, intent(in) :: suppress
        suppress_warnings_after_max = suppress
    end subroutine
    
    subroutine set_scope(self,scope)
        class(exception), intent(inout) :: self
        character(len=*), intent(in) :: scope
        self% scope = scope
    end subroutine set_scope
    
    subroutine set_message(self,message)
        class(exception), intent(inout) :: self
        character(len=*), intent(in) :: message
        self% message = message
    end subroutine set_message
    
    function failure_raise(self,ierr,message) result(raise)
        class(failure), intent(inout) :: self
        integer, intent(in) :: ierr
        character(len=*), intent(in), optional :: message
        logical :: raise
        raise = .FALSE.
        if (ierr == 0) return
        
        raise = .TRUE.
        if (present(message)) self% message = message
        call report_error(self)
    end function failure_raise
    
    function warning_raise(self,ierr,message) result(raise)
        class(warning), intent(inout) :: self
        integer, intent(in) :: ierr
        character(len=*), intent(in), optional :: message
        logical :: raise, raised
        
        raise = .FALSE.
        if (ierr == 0) return
        
        raise = .TRUE.
        warning_count = warning_count + 1
        if (warning_count > max_warnings .and. suppress_warnings_after_max) &
        &   return
        if (present(message)) self% message = message
        call report_error(self)
        if (warning_count == max_warnings .and. suppress_warnings_after_max)  &
        &   write(error_unit,indent4) '(further alerts will be suppressed)'
    end function warning_raise
    
    subroutine set_level(self,level)
        class(alert), intent(inout) :: self
        integer, intent(in) :: level
        self% level = level
    end subroutine set_level
        
    subroutine alert_message(self, msg)
        class(alert), intent(inout) :: self
        character(len=*), intent(in) :: msg
        self% message = msg
        if (verbosity >= self% level) then
            write(error_unit,alert_fmt) trim(self% scope),trim(self% message)
        end if
    end subroutine alert_message

    subroutine set_halt_on_assertion_failure(self,halt)
        class(assertion), intent(inout) :: self
        logical,intent(in) :: halt
        self% halt = halt
    end subroutine set_halt_on_assertion_failure
    
    subroutine assert(self,condition)
        class(assertion), intent(inout) :: self
        logical, intent(in) :: condition
        if (condition) return
        write (error_unit,'(a)') trim(self% scope)//': assertion failed'
        if (self% halt) stop 1
    end subroutine assert
    
    subroutine report_error(test)
        class(exception), intent(inout) :: test
        character(len=7) :: flag
        integer :: verbosity_level

        select type(test)
        type is (failure)
            flag = 'FAILURE'
            verbosity_level = 0
        type is (warning)
            flag = 'Warning'
            verbosity_level = 1
        class default
            flag = 'alert'
            verbosity_level = 1
        end select
        if (verbosity_level <= verbosity) then
            write (error_unit,alert_fmt) flag, trim(test% scope)
            if (len_trim(test% message) > 0)  &
            &   write(error_unit,indent4) trim(test% message)
        end if
    end subroutine report_error
    
end module exceptions_lib
