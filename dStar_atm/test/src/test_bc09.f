program test_bc09
    use constants_def
    use constants_lib
    use nucchem_lib
    use dStar_eos_lib
    use conductivity_lib
    use bc09

    integer :: eos_handle, ierr, i
    real(dp) :: grav, lnTeff, Teff, tau, rho_ph, P_ph
    
    call constants_init('',ierr)
    call nucchem_init('../../data',ierr)
    call dStar_eos_startup('../../data')
    eos_handle = alloc_dStar_eos_handle(ierr)

    grav = 2.43e14_dp
    tau = 2.0*onethird
    rho_ph = -1.0_dp
    
    do i = 20,1,-1
        lnTeff = log10(5.0e5) + (i-1)/19.0
        Teff = 10.0_dp**lnTeff
        call find_photospheric_pressure(Teff,grav,tau,rho_ph,P_ph,eos_handle,ierr)
        if (ierr /= 0) then
            print *,'error: ierr = ',ierr
            rho_ph = -1.0
            cycle
        end if
        write (*,'(3(a,es11.4))') 'rho_ph = ',rho_ph,'    P_ph   = ',P_ph, &
            & '    ln(rho_ph) = ',log(rho_ph)
    end do
    
    call dStar_eos_shutdown
    call nucchem_shutdown

end program test_bc09
