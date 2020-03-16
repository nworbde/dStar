program test_exception
    use, intrinsic:: iso_fortran_env, only: error_unit, output_unit
    use exceptions_lib
    implicit none
    integer :: ierr
    integer :: v
    character(len=message_length) :: message
    type(failure) :: valerr=failure(scope='test_exception')
    type(warning) :: valwarn=warning(scope='test_exception', &
    &   message='while doing verbosity check')
    type(alert) :: two=alert(scope='test_exception',level=2)
    type(alert) :: three=alert(scope='test_exception',level=3)
    type(assertion) :: valcheck=assertion(scope='test_exception')
    real :: val
    
    do v = -1,3
        call set_verbosity(v)
        write(message,'(a,i0)') 'ierr = ',ierr
        write(output_unit,'(/,a,i0)') 'verbosity = ',v
        flush(output_unit)
        ierr = -1
        if (valerr% raised(ierr,message)) then
        end if
        if (valwarn% raised(ierr)) then
        end if
        call two% report('message two')
        call three% report('message three')
        flush(error_unit)
    end do
    
    call set_suppress_warnings_after_max(.TRUE.)
    call set_max_warnings(8)
    call reset_warnings_count
    call valwarn% set_message('')
    
    write(output_unit,'(/,a)') 'testing suppression of repeated warnings.'
    write(output_unit,'(a)') 'maximum number of warnings is 8'
    flush(output_unit)
    do v = 1,10
        write(message,'(a,i0)') 'v = ',v
        if (valwarn% raised(ierr, message)) then
            cycle
        end if
    end do

    write(output_unit,'(/,a)') 'checking assertions'
    flush(output_unit)
    val = 2
    call valcheck% assert(val < 3)
    call valcheck% set_halt_on_assertion_failure(.FALSE.)
    val = 4
    call valcheck% assert(val < 3)
    
end program test_exception

