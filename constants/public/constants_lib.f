module constants_lib
    use constants_def

    implicit none
contains
    
    subroutine constants_init(dstar_dir_init,mesa_dir_init,ierr)
        use const_lib
        use exceptions_lib
        character(len=*), intent(in) :: dstar_dir_init,mesa_dir_init
        integer, intent(out) :: ierr
        type(assertion) :: constants_init_okay = assertion(scope='constants_init')
        ierr = 0
        
        call constants_init_okay% set_message('initializing dstar constants')
        call do_constants_init(dstar_dir_init,ierr)
        call constants_init_okay% assert(ierr==0)
        call constants_init_okay% set_message('initializing mesa constants')
        call const_init(mesa_dir_init,ierr)
        call constants_init_okay% assert(ierr==0)
    end subroutine constants_init

end module constants_lib
