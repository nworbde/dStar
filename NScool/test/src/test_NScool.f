program test_NScool
    use NScool_def
    use NScool_lib
    use NScool_ctrls_io, only: write_controls
    use, intrinsic :: iso_fortran_env, only: output_unit
    
    character(len=*), parameter :: my_dStar_dir = '../../../dStar'
    character(len=*), parameter :: inlist = 'test_inlist'
    
    integer :: ierr, NScool_id
    
    call NScool_init(my_dStar_dir, ierr)
    call check_okay('NScool_init',ierr)
    
    NScool_id = alloc_NScool(ierr)
    call check_okay('NScool_id',ierr)
    
    call NScool_setup(NScool_id,inlist,ierr)
    call check_okay('NScool_setup',ierr)
    
    call write_controls(output_unit,ierr)
    call NScool_shutdown
    
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
end program test_NScool
