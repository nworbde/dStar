program test_bc09
    use constants_def
    use constants_lib
    use nucchem_lib
    use dStar_eos_lib
    use conductivity_lib
    use bc09

    integer :: eos_handle, ierr
    real(dp) :: grav, Teff, tau, Pphoto
    
    call constants_init('',ierr)
    call nucchem_init('../../data',ierr)
    call dStar_eos_startup('../../data')
    eos_handle = alloc_dStar_eos_handle(ierr)

    Teff = 5.0e6_dp
    grav = 2.43e14_dp
    tau = 2.0*onethird
    call find_photospheric_pressure(Teff,grav,tau,Pphoto,eos_handle,ierr)

    call dStar_eos_shutdown
    call nucchem_shutdown

end program test_bc09
