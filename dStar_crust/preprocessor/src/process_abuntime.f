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

    real(dp), parameter :: lgP_increment = 0.005_dp
    character(len=*), parameter ::  &
        & abuntime_filename ='../data/abuntime_lx2_5.data', &
        & abuntime_cache = '../data/cache/abuntime_lx2_5_cache.bin'

    integer :: nz, nion !, ncharged
    real(dp), dimension(:,:), allocatable :: Yion,Yout,Ytot
    character(len=iso_name_length), dimension(:), allocatable :: isos,isonet, &
    &   network
    real(dp) :: T
    integer :: k,ierr, eos_handle, iDelta(1)
    real(dp), dimension(:), allocatable :: lgP,lgPtot
    integer :: nnet, i, j, ntot,nztot

    ierr = 0
    call constants_init('',ierr)
    call check_okay('constants_init',ierr)
    call nucchem_init('../../data/',ierr)
    call check_okay('nucchem_init',ierr)
    
    write(output_unit,'(a)') 'preprocessing abuntime tables'

    write(error_unit,'(/,a)') 'reading abuntime'
    call read_abuntime(abuntime_filename, nz, nion, T, lgP, &
    &   isos, Yion, ierr,lgP_increment)
    call check_okay('read_abuntime',ierr)
    write(error_unit,'(a)') 'done'

    write(error_unit,'(/,a)') 'processing abuntime...'
    call reduce_abuntime(nz,nion,isos,lgP,Yion,nnet,isonet,Yout,ierr)
    call expand_abuntime(nz,nnet,isonet,lgP,Yout,lgP_increment, &
        & nztot,ntot,network,lgPtot,Ytot,ierr)
    write(error_unit,'(a)') 'done'

    write(error_unit,'(/,a)') 'writing cache to '//abuntime_cache//'...'
    call write_abuntime_cache(abuntime_cache,nztot,ntot,network,T, &
    &   lgPtot,Ytot,ierr)
    call check_okay('write_abuntime_cache',ierr)
    write(error_unit,'(a)') 'done'

    deallocate(lgP,isos,Yion,isonet,Yout,lgPtot,Ytot,network)
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
