program test_chem
	use chem_def
	use chem_lib
    use netJina_def, only: iso_name_length

	integer :: j, ierr, ncharged
	integer, dimension(2) :: indx, Z, A, N, charged_ids
	real(dp), dimension(2) :: X, Yion
	character(len=iso_name_length), dimension(2) :: names
	real(dp) :: Xsum
	type(composition_info_type) :: comp
	
	call chem_init('../../data',ierr)
	if (failure('unable to initialize nucchem')) stop
	
	! check the composition average
	Z = [26,28]
	A = 56
    N = A-Z
    
    indx = [(get_nuclide_index_from_ZN(Z(j),N(j)),j=1,2)]
	names = nuclib% name(indx)
    
	X = 0.5
	call compute_composition_moments(2,indx,X,comp,Xsum,  &
	& ncharged, charged_ids, Yion,  &
	& renormalize_mass_fractions = .FALSE., &
	& abunds_are_mass_fractions=.true., exclude_neutrons=.FALSE.)

	write (*,'(/,a,/,18("="))') 'composition moments'
	write (*,'(a10,f8.3)')'<A> = ',comp% A
	write (*,'(a10,f8.3)')'<Z> = ',comp% Z
	write (*,'(a10,f8.3)')'<Z2> = ',comp% Z2
	write (*,'(a10,f8.3)')'Ye = ',comp% Ye
	write (*,'(a10,f8.3)')'<Q> = ',comp% Q
		
	call chem_shutdown
	
	contains
	function failure(str)
		use iso_fortran_env, only : error_unit
		character(len=*), intent(in) :: str
		logical :: failure

        failure = (ierr /= 0)
        if (failure) write (error_unit,'(a)') str
	end function failure
end program test_chem
