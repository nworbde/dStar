program process_abuntime
    use iso_fortran_env, only: error_unit, output_unit
    use constants_def, only: dp
    use constants_lib
    use nucchem_def
    use nucchem_lib
	use superfluid_def
	use superfluid_lib
	use dStar_eos_lib
    use abuntime
    use dStar_crust_mod

    character(len=*), parameter ::  &
        & abuntime_filename ='../data/abuntime_lx2_5.data', &
        & abuntime_cache = '../data/cache/abuntime_lx2_5_cache.bin'

    integer :: nz, nion, ncharged
    real(dp), dimension(:), allocatable :: rho, P, Xneut
    real(dp), dimension(:,:), allocatable :: Yion
    character(len=iso_name_length), dimension(:), allocatable :: isos ! nion
    integer, dimension(:), allocatable :: charged_ids
    type(composition_info_type), dimension(:), allocatable :: ion_info
    real(dp) :: T
    integer :: k,ierr, eos_handle, iDelta(1)
    real(dp), dimension(:), allocatable :: lgP, lgRho, lgEps, delta_lgRho
    real(dp), parameter :: delta = 0.001
    
    ierr = 0
    call constants_init('',ierr)
    call check_okay('constants_init',ierr)
    call nucchem_init('../../data/',ierr)
    call check_okay('nucchem_init',ierr)
	call sf_startup('../../data',ierr)
    call check_okay('sf_startup',ierr)
	call sf_load_gaps('ns','gc','t72',ierr)
    call check_okay('load gaps',ierr)
	
	call dStar_eos_startup('../../data')
	eos_handle = alloc_dStar_eos_handle(ierr)
    call check_okay('allocate eos',ierr)
    
    write(output_unit,'(a)') 'preprocessing abuntime tables'

    write(error_unit,'(/,a)') 'reading abuntime'
    call read_abuntime(abuntime_filename, nz, nion, ncharged, P, rho, T, &
    &   isos, Yion, Xneut, charged_ids, ion_info, ierr)
    call check_okay('read_abuntime',ierr)
    write(error_unit,'(a)') 'done'

    write(error_unit,'(/,a)') 'processing abuntime data...'
    allocate(lgP(nz), lgRho(nz), lgEps(nz), delta_lgRho(nz))
    lgP = log10(P)

    call find_densities(eos_handle, &
    &   lgP,lgRho,lgEps,Yion,ncharged,charged_ids,ion_info,T)

    delta_lgRho = abs(lgRho - log10(rho))
    iDelta = maxloc(delta_lgRho)
    write(error_unit,'(a)') 'done'
    write(error_unit,'(/,a,/,a,f8.5,a,f8.3,", ",f8.3)')  &
    &   'consistency check:','max Delta(lgRho) = ', &
    &   maxval(delta_lgRho),' at lgP, lgRho = ',lgP(iDelta),log10(rho(iDelta))

    write(error_unit,'(/,a)') 'writing cache to '//abuntime_cache//'...'
    call write_abuntime_cache(abuntime_cache,nz,nion,ncharged,isos, &
    &   charged_ids,ion_info,Xneut,T,lgP,lgRho,lgEps,Yion,ierr)
    call check_okay('write_abuntime_cache',ierr)
    write(error_unit,'(a)') 'done'

    deallocate(lgP, lgRho, lgEps, delta_lgRho)
    deallocate(ion_info, charged_ids, isos, Yion, rho, P, Xneut)
    
    call dStar_eos_shutdown
    call sf_shutdown    
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

end program process_abuntime
