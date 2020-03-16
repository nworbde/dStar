module constants_lib
    use constants_def

    implicit none
contains
    
    subroutine constants_init(mesa_dir_init,ierr)
        use const_lib
        use exceptions_lib
        character(len=*), intent(in) :: mesa_dir_init
        integer, intent(out) :: ierr
        type(assertion) :: mesa_const_init = assertion(scope='constants_init')
        ierr = 0
        call const_init(mesa_dir_init,ierr)
        call mesa_const_init% assert(ierr==0)
        call initialize_constants
    end subroutine constants_init

end module constants_lib
