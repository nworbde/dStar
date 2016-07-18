program test_abuntime
    use constants_def, only: dp
    use constants_lib
    use nucchem_def
    use nucchem_lib
	use superfluid_def
	use superfluid_lib
	use dStar_eos_lib
    use abuntime

    character(len=*), parameter ::  &
        & abuntime_filename ='../data/abuntime_lx2_5.data'

    integer :: nz, nion
    real(dp), dimension(:), allocatable :: rho, T
    real(dp), dimension(:,:), allocatable :: Yion, abunds
    character(len=iso_name_length), dimension(:), allocatable :: isos
    integer :: k,ierr
    real(dp), dimension(:), allocatable :: lgP, Xneut
    integer, dimension(:), allocatable :: charged_ids
    integer :: ncharged
    type(composition_info_type), dimension(:), allocatable :: ion_info
    integer :: eos_handle,sf_handle
    
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
    
    print *,'reading abuntime'
    call read_abuntime(abuntime_filename, nz, nion, rho, T, isos, abunds, ierr)
    print *,'done'
    print *,rho(nz),(isos(k),abunds(k,nz),k=1,15)
    
    allocate(Yion(nion,nz),Xneut(nz),charged_ids(nion),ion_info(nz), lgP(nz))
    
    call make_abuntime_crust(eos_handle, rho, T, isos, abunds, lgP, Yion, Xneut, &
    & charged_ids, ncharged, ion_info)

    do k = 1,nz
        write(*,'(2es11.4,3f7.2,es11.4,f7.3)') lgP(k),rho(k),ion_info(k)% Z, ion_info(k)% A, ion_info(k)% Q, Xneut(k), sum(abunds(:,k))
    end do
    
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

end program test_abuntime
