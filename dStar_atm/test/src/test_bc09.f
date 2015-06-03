program test_bc09
    use constants_def
    use constants_lib
    use nucchem_lib
    use dStar_eos_lib
    use conductivity_lib
    use bc09
    use pcy97
    use dStar_atm_def

    integer :: ierr, i
    integer, parameter :: Ntab = 128
    real(dp) :: grav, lgyb, lgy_light, rho_ph, P_ph,kappa
    real(dp) :: lgTeffmin,delta_lgTeff
    real(dp) :: Plight, Pb, lgTb(Ntab), lgTeff_bc(Ntab), lgflux_bc(Ntab)
    real(dp) :: lgTeff_pcy(Ntab), lgflux_pcy(Ntab)
    
    call constants_init('',ierr)
    call nucchem_init('../../data',ierr)
    call dStar_eos_startup('../../data')

    grav = 2.43e14_dp
    lgyb = log10(4.3e13_dp)
    lgy_light = 9.0_dp
    
    Plight = 10.0**lgy_light * grav
    Pb = 10.0**lgyb * grav
    
    lgTeffmin = default_lgTeff_min
    delta_lgTeff = default_lgTeff_max - default_lgTeff_min
    lgTeff_bc = [ (lgTeffmin + real(i-1,dp)*(delta_lgTeff)/real(Ntab-1,dp), &
    &   i = 1, Ntab) ]
    call do_get_bc09_Teff(grav, Plight, Pb, lgTb, lgTeff_bc, lgflux_bc, ierr)
    
    if (ierr /= 0) then
        print *, 'something went wrong: ierr = ', ierr
        stop
    end if

    call do_get_pcy97_Teff(grav, Plight, lgTb, lgTeff_pcy, lgflux_pcy)

    write (*,'(4a10)') 'lg(Te)_bc','lg(Tb)','lg(Te)_pcy','del(Te)'
    write (*,'(4f10.4)') (lgTeff_bc(i),lgTb(i),lgTeff_pcy(i), &
        &   ln10*(lgTeff_bc(i)-lgTeff_pcy(i)),i=1,Ntab)
    
    call dStar_eos_shutdown
    call nucchem_shutdown

end program test_bc09
