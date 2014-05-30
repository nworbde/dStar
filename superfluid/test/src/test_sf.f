program test_sf
	use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
	use superfluid_def
	use superfluid_lib
	
	integer, parameter :: Np = 40
	integer :: i, ierr
	real(dp) :: kF_1, kF_2, k, Tc(max_number_sf_types)
	
	call sf_startup('../../data',ierr)
    if (failure('startup')) stop
	
	call sf_load_gaps('ccdk93','gc','ao85',ierr)
    if (failure('load gaps')) stop
	
	kF_1 = min(sf_tables(neutron_1S0)% kF_min, &
			& sf_tables(neutron_3P2)% kF_min, sf_tables(proton_1S0)% kF_min)
	kF_2 = max(sf_tables(neutron_1S0)% kF_max, &
			&	sf_tables(neutron_3P2)% kF_max, sf_tables(proton_1S0)% kF_max)
	
	write (output_unit,'(a7,tr8,3(a13,tr2))') 'k/fm^-1','Tc(p 1S0)/GK','Tc(n 1S0)/GK','Tc(n 3P2)/GK'
	write (output_unit,'(7("="),tr8,3(13("="),tr2))')
	do i = 1, Np
		k = kF_1 + (kF_2-kF_1)*real(i-1)/real(Np-1)
		call sf_get_results(k,k,Tc)
		write (output_unit,'(f7.3,tr8,3(tr3,f10.6,tr2))') k, Tc*1.0e-9
	end do
	
	call sf_shutdown
	
	contains
        function failure(message)
            character(len=*), intent(in) :: message
            logical :: failure
            
            failure = (ierr /= 0)
            if (failure)  &
            &   write (error_unit,'(a)') 'failure in '//message
        end function failure
	
end program test_sf
