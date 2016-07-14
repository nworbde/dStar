program test_abuntime
    use constants_def, only: dp
    use constants_lib
    use nucchem_def
    use nucchem_lib
    use abuntime

    character(len=*), parameter ::  &
        & abuntime_filename ='../data/abuntime_lx2_5.txt'

    integer :: nz, nion
    real(dp), dimension(:), allocatable :: rho
    real(dp), dimension(:,:), allocatable :: Yion
    character(len=iso_name_length), dimension(:), allocatable :: isos
    integer :: k,ierr
    
    ierr = 0
    call constants_init('',ierr)
    call check_okay('constants_init',ierr)
    
    print *,'reading abuntime'
    call read_abuntime(abuntime_filename, nz, nion, rho, isos, Yion, ierr)
    print *,'done'
    print *,rho(nz),(isos(k),Yion(k,nz),k=1,15)

contains
	subroutine check_okay(msg,ierr)
		use iso_fortran_env, only : error_unit
		character(len=*), intent(in) :: msg
		integer, intent(inout) :: ierr
		if (ierr /= 0) then
			write (error_unit,*) trim(msg)//': ierr = ',ierr
			if (ierr < 0) stop
		end if
	end subroutine check_okay

end program test_abuntime
