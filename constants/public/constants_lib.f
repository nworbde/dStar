module constants_lib
    use constants_def

    implicit none
contains
    
    subroutine constants_init(mesa_dir_init,ierr)
        use const_lib
        character(len=*), intent(in) :: mesa_dir_init
        integer, intent(out) :: ierr
        ierr = 0
        call do_const_init(mesa_dir_init,ierr)
        call initialize_constants
    end subroutine constants_init

end module constants_lib
