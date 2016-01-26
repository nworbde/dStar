program compare_EOS
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

	real(dp), parameter :: lgP_blend_max = 0.35
	real(dp), parameter :: lgP_blend_min = -0.80

	character(len=*), parameter :: datadir='../../data'
	integer, parameter :: Ntrial = 20
	integer :: i, ierr
	real(dp), dimension(Ntrial) :: lgP,lgRho,lgEps,Xprot,dlgRho,dlgEps,dXprot,Xneut
	real(dp), dimension(Ntrial) :: lgPc,lgRhoc,lgEpsc,dlgRhoc,dlgEpsc, fac
    real(dp), dimension(:,:), allocatable :: Yion
    integer, dimension(:), allocatable :: charged_ids
    integer :: ncharged
    type(composition_info_type), dimension(Ntrial) :: ion_info
    integer :: eos_handle, Niso
    real(dp) :: Tref
	
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
	
	lgP = [ (-1.0_dp + 2.0_dp*real(i-1.0,dp)/real(Ntrial-1,dp), i=1,Ntrial) ]

	write (*,'("#",7(a9,tr2),/,"#",7("=========",tr2))') 'lgP','lgRho','dlgRho', &
		& 'lgEps','dlgEps','Xprot','dXprot'

	call do_one('s7a')
	call do_one('s7r')
	call do_one('s7d')
	
    Niso = dStar_crust_get_composition_size()
    allocate(Yion(Niso,Ntrial),charged_ids(Niso))
    
	! convert lgP to cgs
	lgPc = lgP + log10(pressure_n)
    call dStar_crust_get_composition(lgPc,ncharged,charged_ids,Yion,Xneut,ion_info,ierr)
    call check_okay('dStar_crust_get_composition',ierr)
	
    write (*,'("#",a)') 'crust'
	
	do i = 1,Ntrial
        call dStar_crust_get_results(lgPc(i),lgRhoc(i),dlgRhoc(i),lgEpsc(i),dlgEpsc(i),ierr)
        write(*,'(tr1,7(f9.5,tr2))') lgPc(i)-log10(pressure_n), &
        	& lgRhoc(i)-log10(amu*density_n),dlgRhoc(i), &
        	& lgEpsc(i)-log10(mass_n*density_n),dlgEpsc(i),ion_info(i)% Ye,-9.99
    end do
    
	! now do a blended EOS
	where (lgP <= lgP_blend_min)
		lgRho = lgRhoc-log10(amu*density_n)
		lgEps = lgEpsc-log10(mass_n*density_n)
	end where
	where(lgP < lgP_blend_max .and. lgP > lgP_blend_min)
		fac = (lgP-lgP_blend_min)/(lgP_blend_max-lgP_blend_min)
		fac = 0.5_dp*(1.0_dp-cos(pi*fac))
		lgRho = fac*lgRho + (1.0-fac)*(lgRhoc-log10(amu*density_n))
		lgEps = fac*lgEps + (1.0-fac)*(lgEpsc-log10(mass_n*density_n))
	end where

    write (*,'("#",a)') 'blend'
	do i = 1, 20
		write (*,'(tr1,7(f9.5,tr2))') lgP(i),lgRho(i),dlgRho(i), &
			& lgEps(i),dlgEps(i),Xprot(i),dXprot(i)
	end do

	call dStar_core_shutdown
    call dStar_crust_shutdown
    call sf_shutdown
    call nucchem_shutdown

contains
	subroutine do_one(model)
		character(len=*), intent(in) :: model
		integer :: i
		
        write (*,'("#",a)') trim(model)
        
    	call dStar_core_load_table(trim(model),ierr)
    	call check_okay('dStar_core_load_table',ierr)
        
    	do i = 1, 20
    		call dStar_core_get_results(lgP(i),lgRho(i),dlgRho(i), &
    			& lgEps(i),dlgEps(i),Xprot(i),dXprot(i),ierr)
    		call check_okay('dStar_core_get_results',ierr)
    		write (*,'(tr1,7(f9.5,tr2))') lgP(i),lgRho(i),dlgRho(i), &
    			& lgEps(i),dlgEps(i),Xprot(i),dXprot(i)
    	end do
    	call dStar_core_free_table
	end subroutine do_one
	
	subroutine check_okay(msg,ierr)
		use iso_fortran_env, only : error_unit
		character(len=*), intent(in) :: msg
		integer, intent(inout) :: ierr
		if (ierr /= 0) then
			write (error_unit,*) trim(msg)//': ierr = ',ierr
			if (ierr < 0) stop
		end if
	end subroutine check_okay

end program compare_EOS
