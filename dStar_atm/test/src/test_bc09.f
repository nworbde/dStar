program test_bc09
    use constants_def
    use constants_lib
    use nucchem_lib
    use dStar_eos_lib
    use conductivity_lib
    use bc09

    integer :: eos_handle, ierr, i
    real(dp) :: grav, lgTeff, lgyb, lgy_light, lgTb, rho_ph, P_ph,kappa
    
    call constants_init('',ierr)
    call nucchem_init('../../data',ierr)
    call dStar_eos_startup('../../data')
    eos_handle = alloc_dStar_eos_handle(ierr)

    grav = 2.43e14_dp
    lgyb = log10(4.3e13_dp)
    lgy_light = 5.0_dp
    
    write (*,'(5a10)') 'lg(Teff)', 'rho','P/g','kappa','lg(Tb)'
    rho_ph = -1.0_dp
    do i = 30,1,-1
        lgTeff = 5.5_dp + 1.2*(i-1)/29.0
        call do_integrate_bc09_atm(grav,lgyb,lgy_light,lgTeff,lgTb,eos_handle,ierr,rho_ph,P_ph,kappa)
        if (ierr /= 0) then
            print *,'error: ierr = ',ierr
            cycle
        end if
        write (*,'(5f10.4)') lgTeff,rho_ph,P_ph/grav,kappa,lgTb
    end do
    
    call dStar_eos_shutdown
    call nucchem_shutdown

end program test_bc09
