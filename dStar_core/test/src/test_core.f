program test_core
	use constants_def
	use constants_lib
	use nucchem_def
	use nucchem_lib
	use superfluid_lib
	use dStar_eos_lib
	use dStar_core_def
	use dStar_core_lib

	integer, parameter :: Ntrial = 20
	integer :: i, ierr
	real(dp) :: Tref
	real(dp), dimension(Ntrial) :: lgP,lgRho,lgEps,Xprot,dlgRho,dlgEps,dXprot
	integer :: eos_handle
	
	call constants_init('',ierr)
	call check_okay('constants_init',ierr)
	
    call nucchem_init('../../data',ierr)
    call check_okay('nucchem_init',ierr)
        
    call sf_startup('../../data',ierr)
    call check_okay('sf_startup',ierr)
        
    call sf_load_gaps('ns','gc','t72',ierr)
    call check_okay('sf_load_gaps',ierr)
        
    call dStar_eos_startup('../../data')
    call check_okay('dStar_eos_startup',ierr)
	
	call dstar_core_startup('../../data',ierr)
	call check_okay('dstar_core_startup',ierr)

    eos_handle = alloc_dStar_eos_handle(ierr)
    call check_okay('alloc_dStar_eos_handle',ierr)

	Tref = 1.0e8_dp
	call dStar_eos_set_controls(eos_handle,skyrme_parameter_set='s7d')
	call dStar_core_load_table(eos_handle,Tref,ierr)
	call check_okay('dStar_core_load_table',ierr)

	lgP = [(-0.5_dp + 3.5_dp*real(i-1.0,dp)/real(Ntrial-1.0,dp),i=1,Ntrial)]
	lgP = lgP + log10(pressure_n)
	write (*,'(7(a9,tr2),/,7("=========",tr2))') 'lgP','lgRho','dlgRho', &
		& 'lgEps','dlgEps','Xprot','dXprot'
		
    	do i = 1, Ntrial
    		call dStar_core_get_results(lgP(i),lgRho(i),dlgRho(i), &
			&	lgEps(i),dlgEps(i),Xprot(i),dXprot(i),ierr)
    		call check_okay('dStar_core_get_results',ierr)
    		write (*,'(7(f9.5,tr2))') &
    		& lgP(i),lgRho(i),dlgRho(i),lgEps(i),dlgEps(i),Xprot(i),dXprot(i)
    	end do
    	call dStar_core_free_table

	call dStar_core_shutdown
	call dStar_eos_shutdown
	call sf_shutdown
	call nucchem_shutdown
	
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

end program test_core
