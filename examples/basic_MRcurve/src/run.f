program run_MRcurve
	use constants_def
	use constants_lib
    use nucchem_def
    use nucchem_lib
    use superfluid_lib
    use dStar_eos_lib
    use dStar_crust_def
    use dStar_crust_lib
	use dStar_core_def
	use dStar_core_lib
	use dStar_core_tov

	integer :: NMR
	character(len=*), parameter :: datadir='../../data'
	character(len=16) :: model
	integer :: ierr, i, unitno
    integer :: eos_handle
    real(dp) :: Tref, lgPstart, lgPend, lnP_rez
	real(dp), dimension(:), pointer :: y
	real(dp), dimension(:), allocatable :: M, R, P_c, rho_c, eps_c
	real(dp) :: lgP, lgRho, lgEps, lgPmin, lgPmax
	
	namelist /MRcontrols/ model, NMR, lgPmin, lgPmax
	
	open(newunit=unitno,file='inlist',status='old',action='read') 
	read(unitno,nml=MRcontrols,iostat=ierr)
	close(unitno)
	
	if (ierr /= 0) then
		write(*,*) 'unable to read controls'
		stop
	end if
	
	call constants_init('',ierr)
	call check_okay('constants_init',ierr)

	call dstar_core_startup(datadir,ierr)
	call check_okay('dstar_core_startup',ierr)
    
    call nucchem_init(datadir,ierr)
    call check_okay('nucchem_init',ierr)
	
    call sf_startup(datadir,ierr)
    call check_okay('sf_startup',ierr)
	
    call sf_load_gaps('ns','gc','t72',ierr)
    call check_okay('sf_load_gaps',ierr)
	
    call dStar_eos_startup(datadir)
    call check_okay('dStar_eos_startup',ierr)
	
    eos_handle = alloc_dStar_eos_handle(ierr)
    call check_okay('alloc_dStar_eos_handle',ierr)
    
    ! switch off the warnings about quantum effects
    call dStar_eos_set_controls(eos_handle,suppress_warnings=.TRUE.)
    
    call dStar_crust_startup(datadir,ierr)
    call check_okay('dStar_atm_startup',ierr)
    
    Tref = 1.0d8
    call dStar_crust_load_table('hz90',eos_handle, Tref,ierr)
    call check_okay('dStar_crust_load_table',ierr)
	
	call dStar_core_load_table(model,ierr)
	call check_okay('dStar_core_load_table',ierr)
	
	lgPend = -9.0
	lnP_rez = 0.2
	allocate(y(num_core_variables))
	allocate(M(NMR),R(NMR),P_c(NMR),rho_c(NMR),eps_c(NMR))
	
	do i = 1,NMR
		lgPstart = lgPmin + (lgPmax-lgPmin)*real(i-1,dp)/real(NMR-1,dp)
		call core_integrate(lgPstart, lgPend, lnP_rez, y, ierr)
		call check_okay('core_integrate',ierr)
		M(i) = y(core_mass)
		R(i) = y(core_radius)*length_g*1.0e-5
		P_c(i) = core_structure% pressure(1) * pressure_g/pressure_n
		lgP = log10(P_c(i))
		call core_get_EOS(lgP,lgRho,lgEps,ierr)
		rho_c(i) = 10.0_dp**lgRho
		eps_c(i) = 10.0_dp**lgEps
		call core_write_crust('LOGS',model)
	end do

	open(newunit=unitno,file='MRcurve_'//trim(model),action='write')
	write(unitno,'(5a10)') 'M (Msun)','R (km)','P_c','n/n0','eps_c'
	do i = 1, NMR
		write (unitno,'(5f10.4)') M(i), R(i), P_c(i), rho_c(i)/0.16, eps_c(i)
	end do
	close(unitno)
	
	deallocate(y)
	deallocate(M,R,P_c,rho_c,eps_c)

    call dStar_crust_shutdown
	call dStar_eos_shutdown
    call sf_shutdown
    call nucchem_shutdown
	call dStar_core_shutdown
	
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

end program run_MRcurve
