program test_nucchem
    use exceptions_lib
    use constants_lib
	use nucchem_def
	use nucchem_lib

	integer :: j, ierr, ncharged
	integer, dimension(2) :: indx, Z, A, N, charged_ids
	real(dp), dimension(2) :: X, Yion
	character(len=iso_name_length), dimension(2) :: names
	real(dp) :: Xsum
	type(composition_info_type) :: comp
    type(assertion) :: initialization=assertion(scope='main', &
    &   message='nucchem is uninitialized')
	
    call constants_init('../..','',ierr)
    call initialization% assert(ierr==0)
    call nucchem_init(ierr)
    call initialization% assert(ierr==0)
	
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
		
	call nucchem_shutdown
end program test_nucchem
