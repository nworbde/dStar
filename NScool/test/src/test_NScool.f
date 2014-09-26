program test_NScool
    use NScool_def
    use NScool_lib
    use create_model
    use NScool_ctrls_io, only: write_controls
    use, intrinsic :: iso_fortran_env, only: output_unit
    
    character(len=*), parameter :: my_dStar_dir = '../../../dStar'
    character(len=*), parameter :: inlist = 'test_inlist'
    type(NScool_info), pointer :: s
    integer :: ierr, NScool_id,i
    
    call NScool_init(my_dStar_dir, ierr)
    call check_okay('NScool_init',ierr)
    
    NScool_id = alloc_NScool(ierr)
    call check_okay('NScool_id',ierr)
    
    call NScool_setup(NScool_id,inlist,ierr)
    call check_okay('NScool_setup',ierr)
    ierr = 0
    call write_controls(output_unit,ierr)
    
    call get_NScool_info_ptr(NScool_id,s,ierr)
    call do_setup_crust_zones(s, ierr)
    
    call tov_write_crust
    do i = s% nz, 1, -1
        write (output_unit,'(es15.8,tr2,f14.10,tr1,es15.8,f14.10)') s% P_bar(i), s% ePhi_bar(i), s% m(i), s% eLambda_bar(i)
    end do
    
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
