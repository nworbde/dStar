module abuntime
    ! reader for abuntime files
    use const_def, only: dp
    
    real(dp), parameter :: gravity = 1.85052e14_dp
    real(dp), parameter :: mdot = 8.8e4_dp*0.3_dp
    real(dp), parameter :: default_min_P_increment = 0.0_dp
    
contains

    subroutine read_abuntime(abuntime_filename, nz, nion, ncharged, P, rho, T, &
    &   isos, Yion, Xneut, charged_ids, ion_info, ierr, min_P_increment)
        
        use iso_fortran_env, only: iostat_end, error_unit
        use const_def
        use nucchem_def
        use nucchem_lib
        
        character(len=*), intent(in) :: abuntime_filename
        integer, intent(out) :: nz  ! number of zones
        integer, intent(out) :: nion    ! number of species
        integer, intent(out) :: ncharged
        real(dp), intent(out), dimension(:), allocatable :: P, rho ! nz
        real(dp), intent(out) :: T ! reference temperature
        character(len=iso_name_length), dimension(:), allocatable :: isos ! nion
        real(dp), intent(out), dimension(:,:), allocatable :: Yion   ! (nion,nz)
        real(dp), dimension(:), intent(out), allocatable :: Xneut
        integer, intent(out), dimension(:), allocatable :: charged_ids
        type(composition_info_type), dimension(:), intent(out), allocatable :: ion_info
        integer, intent(out) :: ierr
        real(dp), intent(in), optional :: min_P_increment
 
        integer, dimension(max_nnuclib) :: network_indcs
        integer, dimension(:), allocatable :: indcs
        real(dp) :: Xsum
        integer :: i, n_indx
 
        integer, parameter :: default_chunk_size = 4096
        integer :: unitno, ios
        integer :: k,nrejected
        real(dp), dimension(:), allocatable :: abunds
        real(dp) :: time,temp,EFe,EFn,last_P,delta_P
        
        if (present(min_P_increment)) then
            delta_P = 1.0+min_P_increment
        else
            delta_P = 1.0+default_min_P_increment
        end if
        
        open(newunit=unitno,file=trim(abuntime_filename), &
        &   action='read',status='old', iostat=ierr)
        if (failure('opening'//trim(abuntime_filename),ierr)) return
        
        write(error_unit,'(a)') 'reading '//abuntime_filename
        read(unitno,'(1x,i5)') nion
        write(error_unit,'(a,i5)') 'nion = ',nion
        allocate(isos(nion),abunds(nion),charged_ids(nion))
                
        allocate(P(default_chunk_size), &
        & rho(default_chunk_size),Yion(nion,default_chunk_size), &
        & ion_info(default_chunk_size), Xneut(default_chunk_size), stat=ierr)
        if (failure('allocating P, rho, Yion, ion_info, Xneut',ierr)) return

        nz = 1
        last_P = 0.0
        nrejected = 0
        do
            ! size check
            if (nz > size(rho)) then
                call realloc_abuntime_arrays(2*size(P),ierr)
                if (failure('reallocating abuntime arrays',ierr)) return
            end if
            read(unitno,'(1x,1e18.10,2e11.3,2X,2E18.10)',iostat=ios) &
                 & time,T,rho(nz),EFe,EFn
            if (ios == iostat_end) then
                nz = nz-1
                exit
            else if (ios /= 0) then
                write(error_unit,*) 'abnormal return: ierr = ', ierr
                return
            end if
            P(nz) = gravity*mdot*time       
            read(unitno,'(1x,5(a5,1pe10.3))') (isos(k),abunds(k),k=1,nion)
            read(unitno,*)
            
            if (nz == 1) then
                ! set the network pointers
                allocate(indcs(nion))
                indcs = [(get_nuclide_index(adjustl(isos(i))),i=1,nion)]

                do i = 1, nion
                    if (indcs(i) == -1) then
                        print *,isos(i),': index not found'
                    end if
                end do

                ! find the neutron indx
                network_indcs = 0
                network_indcs(indcs) = [(i,i=1,nion)]
                n_indx = network_indcs(get_nuclide_index('n'))
            end if
            
            call compute_composition_moments(nion,indcs,abunds(:), &
            &   ion_info(nz),Xsum,ncharged,charged_ids,Yion(:,nz), &
            &   exclude_neutrons=.TRUE., abunds_are_mass_fractions=.FALSE.)
                        
            ! only increment nz if P has increased
            if (P(nz) > last_P*delta_P) then
                last_P = P(nz)
                nz = nz+1
            else
                nrejected = nrejected+1
            end if
            if (modulo(nz,1000) == 0)  &
            & write (error_unit,'(a)',advance='no') '.'
        end do
        write(error_unit,'(/,a,i0,a)') 'got ',nz,' zones'
        write(error_unit,'(a,i0,a)') 'rejected ',nrejected,' zones: dP = 0'
        close(unitno)
        
        call realloc_abuntime_arrays(nz,ierr)
        
        ! scale temperature to Kelvin
        T = T* 1.0e9_dp
        
    contains
        subroutine realloc_abuntime_arrays(newsize,ierr)
            integer, intent(in) :: newsize
            integer, intent(out) :: ierr
            integer :: currentsize, oldsize
            real(dp), dimension(:), allocatable :: tmp_rho, tmp_P, tmp_Xneut
            real(dp), dimension(:,:), allocatable :: tmp_Yion
            type(composition_info_type), dimension(:), allocatable :: tmp_ion_info
            
            oldsize = size(rho)
            allocate(tmp_rho(newsize),tmp_Yion(nion,newsize),tmp_P(newsize), &
            &   tmp_ion_info(newsize), tmp_Xneut(newsize), stat=ierr)
            if (ierr /= 0) return
            currentsize = min(oldsize,newsize)
            tmp_rho(1:currentsize) = rho(1:currentsize)
            tmp_P(1:currentsize) = P(1:currentsize)
            tmp_Yion(:,1:currentsize) = Yion(:,1:currentsize)
            tmp_ion_info(1:currentsize) = ion_info(1:currentsize)
            tmp_Xneut(1:currentsize) = Xneut(1:currentsize)
            deallocate(rho,P,Yion,ion_info,Xneut)
            allocate(rho(newsize),P(newsize),Yion(nion,newsize), &
            &   ion_info(newsize), Xneut(newsize), stat=ierr)
            if (ierr /= 0) return
            rho(1:newsize) = tmp_rho(1:newsize)
            P(1:newsize) = tmp_P(1:newsize)
            Yion(:,1:newsize) = tmp_Yion(:,1:newsize)
            ion_info(1:newsize) = tmp_ion_info(1:newsize)
            Xneut(1:newsize) = tmp_Xneut(1:newsize)
            deallocate(tmp_rho,tmp_P,tmp_Yion,tmp_ion_info,tmp_Xneut)
        end subroutine realloc_abuntime_arrays
    end subroutine read_abuntime
        
    subroutine read_abuntime_cache(cache_filename,nz,nion,ncharged,isos, &
        &   charged_ids,ion_info,T,lgP,lgRho,lgEps,Yion,ierr)
        use nucchem_def, only: iso_name_length, composition_info_type

        character(len=*), intent(in) :: cache_filename
        integer, intent(out) :: nz,nion,ncharged
        character(len=iso_name_length), intent(out), dimension(:), allocatable :: isos
        type(composition_info_type), dimension(:), allocatable, intent(out):: &
            &   ion_info
        real(dp), intent(out) :: T
        real(dp), dimension(:), intent(out), allocatable :: lgP, lgRho, lgEps
        integer, intent(out), dimension(:), allocatable :: charged_ids
        real(dp), dimension(:,:), intent(out), allocatable :: Yion
        integer, intent(out) :: ierr
        integer :: unitno
        
        open(newunit=unitno,file=trim(cache_filename), &
        &   action='read',status='old',form='unformatted', iostat=ierr)
        if (failure('opening'//trim(cache_filename),ierr)) return
                
        read(unitno) nz
        read(unitno) nion
        read(unitno) ncharged
        allocate(isos(nion),ion_info(nz), charged_ids(nion), &
        &   lgP(nz),lgRho(nz),lgEps(nz),Yion(nion,nz),stat=ierr)
        if (failure('allocating abuntime tables',ierr)) then
            close(unitno)
            return
        end if
        read(unitno) isos
        read(unitno) charged_ids
        read(unitno) ion_info
        read(unitno) T
        read(unitno) lgP
        read(unitno) lgRho
        read(unitno) lgEps
        read(unitno) Yion
        close(unitno)
              
    end subroutine read_abuntime_cache

    subroutine write_abuntime_cache(cache_filename,nz,nion,ncharged,isos, &
    &   charged_ids,ion_info,T,lgP,lgRho,lgEps,Yion,ierr)

        use nucchem_def, only: iso_name_length, composition_info_type
        character(len=*), intent(in) :: cache_filename
        integer, intent(in) :: nz, nion, ncharged
        character(len=iso_name_length), intent(in), dimension(:) :: isos
        integer, intent(in), dimension(:) :: charged_ids
        type(composition_info_type), intent(in), dimension(:) :: ion_info
        real(dp), intent(in) :: T
        real(dp), intent(in), dimension(:) :: lgP, lgRho, lgEps
        real(dp), intent(in), dimension(:,:) :: Yion
        integer, intent(out) :: ierr
        integer :: unitno
        
        ierr = 0
        open(newunit=unitno, file=trim(cache_filename),action='write', &
        &   form='unformatted',iostat=ierr)
        if (failure('opening '//trim(cache_filename),ierr)) return

        write(unitno) nz
        write(unitno) nion
        write(unitno) ncharged
        write(unitno) isos
        write(unitno) charged_ids
        write(unitno) ion_info
        write(unitno) T
        write(unitno) lgP
        write(unitno) lgRho
        write(unitno) lgEps
        write(unitno) Yion
        close(unitno)
    end subroutine write_abuntime_cache
        
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
