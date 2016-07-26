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

    character(len=*), parameter ::  &
        & abuntime_filename ='../data/abuntime_lx2_5.data', &
        & abuntime_cache = '../data/cache/abuntime_lx2_5_cache.bin'

    integer :: nz, nion, ncharged
    real(dp), dimension(:), allocatable :: rho, P, Xneut, deltaRho
    real(dp), dimension(:,:), allocatable :: Yion
    character(len=iso_name_length), dimension(:), allocatable :: isos ! nion
    integer, dimension(:), allocatable :: charged_ids
    type(composition_info_type), dimension(:), allocatable :: ion_info
    real(dp) :: T
    integer :: k,ierr, eos_handle
    
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
    
    write(output_unit,*) 'preprocessing abuntime tables'

    write(error_unit,*) 'reading abuntime'
    call read_abuntime(abuntime_filename, nz, nion, ncharged, P, rho, T, &
    &   isos, Yion, Xneut, charged_ids, ion_info, ierr)
    call check_okay('read_abuntime',ierr)
    write(error_unit,*) 'done'

    print *,'ncharged = ',ncharged
    print *,'size charged_ids = ',size(charged_ids)

    write(error_unit,*) 'writing cache'
    call write_abuntime_cache(abuntime_cache,nz,nion,ncharged,isos, &
    &   charged_ids,ion_info,P,Yion,ierr)    
    call check_okay('write_abuntime_cache',ierr)
    write(error_unit,*) 'done'

    allocate(deltaRho(nz))
    call check_abuntime(eos_handle, rho, T, isos, P,  &
    &   Yion, Xneut, charged_ids, ncharged, ion_info, deltaRho)
    
    do k = 1, nz
        write(output_unit,'(2(es11.4,tr2),2(f3.0,tr1),f6.3,tr1,f5.2)')  &
        &   P(k), rho(k), ion_info(k)% Z, ion_info(k)% A,  &
        &   ion_info(k)% Q,Xneut(k)
    end do
    
	call nucchem_shutdown
    
    write(error_unit,'(a,es11.4)') 'max delta Rho = ',maxval(deltaRho)
    write(error_unit,'(a,es11.4)') 'at P = ',P(maxloc(deltaRho))
    
    deallocate(deltaRho,ion_info,charged_ids,Xneut,Yion,P,isos,rho)

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
