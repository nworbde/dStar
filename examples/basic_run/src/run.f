program run_dStar
    use NScool_def
    use NScool_lib
    
    character(len=*), parameter :: my_dStar_dir = '/path/to/local/dStar'
    character(len=*), parameter :: inlist = 'inlist'
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
    ! the parameters are in the inlist
    call NScool_evolve_model(NScool_id,ierr)    

    ! now start the cooling
    ! to do this, we need to reset some of our paramters
    
    ! we first get a pointer to the star data structure
    call get_NScool_info_ptr(NScool_id,s,ierr)

    ! Now set Mdot = 0
    s% Mdot = 0.0_dp
    ! set the starting model profile to be one larger than the current model
    s% starting_number_for_profile = s% model + 1
    ! We can reset the start time to zero for convenience
    s% start_time = 0.0
    ! and we'll have a 1000 day quiescent period
    s% maximum_end_time = 8.64d8
    
    ! now evolove our cooling model
    call NScool_evolve_model(NScool_id,ierr)
    
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
end program run_dStar
