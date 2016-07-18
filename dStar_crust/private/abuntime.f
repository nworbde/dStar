module abuntime
    ! reader for abuntime files

contains
    
    subroutine make_abuntime_crust(eos_handle, rho, T, isos, abunds, lgP,  &
    &   Yion, Xneut, charged_ids, ncharged, ion_info)
        ! call this after read_abuntime
        use const_def
        use nucchem_def
        use nucchem_lib
    	use superfluid_def
    	use superfluid_lib
    	use dStar_eos_lib
        
        integer, intent(in) :: eos_handle
        real(dp), dimension(:), intent(in) :: rho
        real(dp), dimension(:), intent(in) :: T
        character(len=iso_name_length), dimension(:), intent(in) :: isos
        real(dp), dimension(:,:), intent(in) :: abunds
        real(dp), dimension(:), intent(out) :: lgP
        real(dp), dimension(:,:), intent(out) :: Yion
        real(dp), dimension(:), intent(out) :: Xneut
        integer, intent(out), dimension(:) :: charged_ids
        integer, intent(out) :: ncharged
        type(composition_info_type), dimension(:), intent(out) :: ion_info
        !
        integer, dimension(max_nnuclib) :: network_indcs
        integer, dimension(:), allocatable :: indcs
        real(dp) :: Xsum
        integer :: nz, nion
        integer :: i, n_indx
		real(dp) :: k,Tcs(max_number_sf_types),chi,kFn,kFp
		integer :: phase
    	type(crust_eos_component), dimension(num_crust_eos_components) :: &
    	& eos_components
    	real(dp), dimension(num_dStar_eos_results) :: res
        
        nz = size(rho)
        nion = size(isos)
        
        allocate(indcs(nion))
        
        ! set the network pointers
        indcs = [(get_nuclide_index(adjustl(isos(i))),i=1,nion)]
        
        do i = 1, nion
            if (indcs(i) == -1) then
                print *,isos(i),': index not found'
            end if
        end do

        ! and the reverse lookup
        network_indcs = 0
        network_indcs(indcs) = [(i,i=1,nion)]
        n_indx = network_indcs(get_nuclide_index('n'))
        
        do i = 1, nz
            call compute_composition_moments(nion,indcs,abunds(:,i), &
            &   ion_info(i),Xsum,ncharged,charged_ids,Yion(:,i), &
            &   exclude_neutrons=.TRUE., abunds_are_mass_fractions=.FALSE.)
            
			chi = nuclear_volume_fraction(rho(i),ion_info(i), &
			&   default_nuclear_radius)
			kFn = neutron_wavenumber(rho(i),ion_info(i),chi)
			kFp = 0.0_dp
			call sf_get_results(kFp,kFn,Tcs)
			call eval_crust_eos(eos_handle,rho(i),T(i),ion_info(i), &
			&   ncharged,charged_ids, Yion(:,i), Tcs, res, phase, chi, &
			&   eos_components)
            
            lgP(i) = res(i_lnP)/ln10
        end do
        Xneut = abunds(n_indx,:)
        
    end subroutine make_abuntime_crust
    
    subroutine read_abuntime(abuntime_filename, nz, nion, rho, T,  &
    &   isos, Yion, ierr)
        
        use iso_fortran_env, only: iostat_end, error_unit
        use const_def
        use nucchem_def, only : iso_name_length
        
        character(len=*), intent(in) :: abuntime_filename
        integer, intent(out) :: nz  ! number of zones
        integer, intent(out) :: nion    ! number of species
        real(dp), intent(out), dimension(:), allocatable :: rho  ! nz
        real(dp), intent(out), dimension(:), allocatable :: T    ! nz
        character(len=iso_name_length), dimension(:), allocatable :: isos ! nion
        real(dp), intent(out), dimension(:,:), allocatable :: Yion   ! (nion,nz)
        integer, intent(out) :: ierr
 
        integer, parameter :: default_chunk_size = 4096
        integer :: unitno, ios
        integer :: k
        real(dp) :: time,temp,EFe,EFn,last_rho
        
        open(newunit=unitno,file=trim(abuntime_filename), &
        &   action='read',status='old', iostat=ierr)
        if (failure('opening'//trim(abuntime_filename),ierr)) return
        
        print *,'reading '//abuntime_filename
        read(unitno,'(1x,i5)') nion
        print *,'nion = ',nion
        allocate(isos(nion))
        
        allocate(rho(default_chunk_size),Yion(nion,default_chunk_size), &
            & T(default_chunk_size), stat=ierr)
        if (failure('allocating rho,Yion',ierr)) return

        nz = 1
        last_rho = 0.0
        do
            ! size check
            if (nz > size(rho)) then
                call realloc_abuntime_arrays(2*size(rho),ierr)
                if (failure('reallocating abuntime arrays',ierr)) return
            end if
            read(unitno,'(1x,1e18.10,2e11.3,2X,2E18.10)',iostat=ios) &
                 & time,T(nz),rho(nz),EFe,EFn
            if (ios == iostat_end) then
                nz = nz-1
                exit
            else if (ios /= 0) then
                print *,'abnormal return: ierr = ', ierr
                return
            end if
            ! only accept if rho has incremented
            read(unitno,'(1x,5(a5,1pe10.3))') (isos(k),Yion(k,nz),k=1,nion)
            read(unitno,*)
            ! only increment nz if rho has increased
            if (rho(nz) > last_rho) then
                last_rho = rho(nz)
                nz = nz+1
            end if
            if (modulo(nz,100) == 0)  &
            & write (error_unit,'(a)',advance='no') '.'
        end do
        print *,'got ',nz,' zones'
        close(unitno)
        
        call realloc_abuntime_arrays(nz,ierr)
        
    contains
        subroutine realloc_abuntime_arrays(newsize,ierr)
            integer, intent(in) :: newsize
            integer, intent(out) :: ierr
            integer :: currentsize, oldsize
            real(dp), dimension(:), allocatable :: tmp_rho, tmp_T
            real(dp), dimension(:,:), allocatable :: tmp_Yion
            
            oldsize = size(rho)
            allocate(tmp_rho(newsize),tmp_T(newsize),tmp_Yion(nion,newsize), &
            &   stat=ierr)
            if (ierr /= 0) return
            currentsize = min(oldsize,newsize)
            tmp_rho(1:currentsize) = rho(1:currentsize)
            tmp_T(1:currentsize) = T(1:currentsize)
            tmp_Yion(:,1:currentsize) = Yion(:,1:currentsize)
            deallocate(rho,T,Yion)
            allocate(rho(newsize),T(newsize),Yion(nion,newsize),stat=ierr)
            if (ierr /= 0) return
            rho(1:newsize) = tmp_rho(1:newsize)
            T(1:newsize) = tmp_T(1:newsize)
            Yion(:,1:newsize) = tmp_Yion(:,1:newsize)
            deallocate(tmp_rho,tmp_T,tmp_Yion)
        end subroutine realloc_abuntime_arrays
    end subroutine read_abuntime

    function failure(msg,ierr)
        use, intrinsic :: iso_fortran_env, only: error_unit
        character(len=*), intent(in) :: msg
        integer, intent(in) :: ierr
        logical :: failure
        
        failure = (ierr /= 0)
        if (failure) then
            write(error_unit,*) trim(msg),': ierr = ',ierr
        end if
    end function failure
    
end module abuntime
