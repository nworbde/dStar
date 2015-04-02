program test_NScool
    use iso_fortran_env, only : output_unit
    use NScool_def
    use NScool_lib

    character(len=*), parameter :: my_dStar_dir = '../../'
    character(len=*), parameter :: inlist = 'test_inlist'
    type(NScool_info), pointer :: s
    integer :: ierr, NScool_id
    
    call NScool_init(my_dStar_dir, ierr)
    call check_okay('NScool_init',ierr)
    
    NScool_id = alloc_NScool(ierr)
    call check_okay('NScool_id',ierr)
    
    call NScool_setup(NScool_id,inlist,ierr)
    call check_okay('NScool_setup',ierr)
    ierr = 0

    call NScool_create_model(NScool_id,ierr)
    
    ! run with accretion on; heat the crust
    call NScool_evolve_model(NScool_id,ierr)    
    call get_NScool_info_ptr(NScool_id,s,ierr)

    ! now start the cooling
    s% Mdot = 0.0_dp
    s% starting_number_for_profile = s% model + 1
    s% start_time = s% tsec
    ! specify new maximum_end_time for length of cooling 
    s% maximum_end_time = 1.0d3*(s% tsec)
    
    call NScool_evolve_model(NScool_id,ierr)
    
    call NScool_shutdown
    
    write (output_unit,'(/,/,a,i3)') 'test_NScool exited with ierr = ',ierr
    
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
