program print_bc09_atm
    use math_lib
    use exceptions_lib
    use constants_lib
    use constants_lib
    use nucchem_lib
    use dStar_eos_lib
    use conductivity_lib
	use dStar_atm_def
	use dStar_atm_lib
    use bc09
	
    character(len=*), parameter :: datadir = '../../data'
    real(dp), parameter :: grav = 2.43e14_dp, Pb = 1.0e27_dp, Plight = 1.0e8_dp*grav, lgTeff = 5.3_dp
    logical, parameter :: print_atm = .TRUE.
	integer :: ierr,i
    integer :: eos_handle, cond_handle
    real(dp) :: lgyb, lgy_light, rho_ph, lgTb
    type(failure) :: integration_failure=failure(scope='print_bc09_atm', &
    &   message='while integrating over atmosphere')
    character(len=128) :: msg
    type(alert) :: status=alert(scope='do_get_bc09_Teff')

    
    call math_init()
    call constants_init('../..','',ierr)
    call check_okay('constants_init',ierr)
    call nucchem_init(ierr)
    call check_okay('nucchem_init',ierr)
    call dStar_eos_startup(ierr)
    call check_okay('eos_startup',ierr)
    call conductivity_startup(ierr)
    call check_okay('conductivity_startup',ierr)
    call dStar_atm_startup(ierr)
    call check_okay('dStar_atm_startup',ierr)
    
    ! constants, nucchem, and dStar_eos must be initialized
    eos_handle = alloc_dStar_eos_handle(ierr)
    cond_handle = alloc_conductivity_handle(ierr)
    call conductivity_set_controls(cond_handle, &
    &   include_neutrons=.FALSE., &
    &   include_superfluid_phonons=.FALSE., &
    &   lgrho_rad_off=10.0_dp, lgrho_rad_on=9.0_dp)
    
    ! set the integration domain
    lgyb = log10(Pb/grav)
    lgy_light = log10(Plight/grav)
    
    rho_ph = -1.0_dp

    ! print header
    write(*,'(4(a11,tr1),2(a4,tr1))') 'density','temperature','pressure','opacity','Z','A'
    write(*,'(4(a11,tr1)//58("="))') 'g cm^-3','K','dyn cm^-2','cm^2 g^-1'
    
    call do_integrate_bc09_atm( &
    &   grav,lgyb,lgy_light,lgTeff,lgTb,rho_ph, &
    &   eos_handle,cond_handle,ierr,print_atm)
    if (integration_failure% raised(ierr)) then
        write(*,'(a,i3)') 'error: ierr = ',ierr
    end if

    call free_conductivity_handle(cond_handle)
    call free_dStar_eos_handle(eos_handle)
    
    call conductivity_shutdown
    call dStar_eos_shutdown
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

end program print_bc09_atm

