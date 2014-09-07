program test_crust
    use constants_def, only: dp
    use constants_lib
    use nucchem_def
    use nucchem_lib
    use superfluid_lib
    use dStar_eos_lib
    use hz90
    
    integer, parameter :: Ntab = 500
    integer :: i, ierr
    real(dp), dimension(Ntab) :: lgP, lgrho
    real(dp), dimension(HZ90_number,Ntab) :: Yion
    integer, dimension(HZ90_number) :: charged_ids
    integer :: ncharged
    real(dp), dimension(Ntab) :: Xneut
    type(composition_info_type), dimension(Ntab) :: ion_info
    integer :: eos_handle
    
    call constants_init('',ierr)
    call nucchem_init('../../data',ierr)
	call sf_startup('../../data',ierr)
	call sf_load_gaps('ns','gc','t72',ierr)
	
	call dStar_eos_startup('../../data')
	eos_handle = alloc_dStar_eos_handle(ierr)
    
    if (ierr /= 0) stop
    
    lgP = [(26.5+5.0*real(i-1,dp)/real(Ntab-1,dp),i=1,Ntab)]
    call do_make_crust(lgP,Yion,Xneut,charged_ids,ncharged,ion_info)
    
    call find_densities(eos_handle,lgP,lgRho,Yion,ncharged,charged_ids,ion_info)
    write(*,'(21a9,/)') 'lg P','lg rho',HZ90_network(:)
    do i = 1,Ntab
        write(*,'(21f9.5)') lgP(i),lgRho(i),Xneut(i),Yion(1:ncharged,i)
        if (mod(i,20) == 0) then
            write(*,'(21a9,/)') 'lg P','lg rho',HZ90_network(:)
        end if  
    end do
    call sf_shutdown
    call nucchem_shutdown

end program test_crust
