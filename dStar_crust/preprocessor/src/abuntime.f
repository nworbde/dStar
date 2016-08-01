module abuntime
    ! reader for abuntime files
    use const_def, only: dp, pi
    use utils_def
    use utils_lib
    use nucchem_def, only: iso_name_length
    
    real(dp), parameter :: gravity = 1.85052e14_dp
    real(dp), parameter :: mdot = 8.8e4_dp*0.3_dp
    real(dp), parameter :: default_lgP_increment = 0.0_dp
    real(dp), parameter :: abundance_threshold = 0.01_dp
    
    type(integer_dict), pointer, save :: isodir=>null()
    
contains

    subroutine read_abuntime(filename, nz, nion, T, lgP, isos, Y, ierr, lgP_increment)
        use iso_fortran_env, only: iostat_end, error_unit
        use const_def
        use nucchem_def

        character(len=*), intent(in) :: filename
        integer, intent(out) :: nz  ! number of zones
        integer, intent(out) :: nion    ! number of species
        real(dp), intent(out) :: T
        real(dp), intent(out), dimension(:), allocatable :: lgP ! nz
        character(len=iso_name_length), dimension(:), allocatable :: isos ! nion
        real(dp), intent(out), dimension(:,:), allocatable :: Y   ! (nion,nz)
        integer, intent(out) :: ierr
        real(dp), intent(in), optional :: lgP_increment

        integer, parameter :: default_chunk_size = 4096
        integer :: unitno, ios
        integer :: k,nrejected,nloop
        real(dp) :: time,rho,EFe,EFn,last_lgP,delta_lgP

        delta_lgP = default_lgP_increment
        if (present(lgP_increment)) delta_lgP = lgP_increment

        open(newunit=unitno,file=trim(filename), &
        &   action='read',status='old', iostat=ierr)
        if (failure('opening'//trim(filename),ierr)) return

        read(unitno,'(1x,i5)') nion
        write(error_unit,'(a,i5)') 'nion = ',nion
        allocate(isos(nion))

        allocate(lgP(default_chunk_size), Y(nion,default_chunk_size), &
        &   stat=ierr)
        if (failure('allocating P, rho, Yion, ion_info, Xneut',ierr)) return

        nz = 1
        nloop = 1
        last_lgP = -100.0
        nrejected = 0
        do
            ! size check
            if (nz > size(lgP)) then
                call realloc_abuntime_arrays(2*size(lgP),ierr)
                if (failure('reallocating abuntime arrays',ierr)) return
            end if
            read(unitno,'(1x,1e18.10,2e11.3,2X,2E18.10)',iostat=ios) &
            & time,T,rho,EFe,EFn
            if (ios == iostat_end) then
                nz = nz-1
                exit
            else if (ios /= 0) then
                write(error_unit,*) 'abnormal return: ierr = ', ierr
                return
            end if
            lgP(nz) = log10(gravity*mdot*time)
            read(unitno,'(1x,5(a5,1pe10.3))') (isos(k),Y(k,nz),k=1,nion)
            read(unitno,*)

            ! only increment nz if lgP has increased
            if (lgP(nz) > last_lgP+delta_lgP) then
                last_lgP = lgP(nz)
                nz = nz+1
            else
                nrejected = nrejected+1
            end if
            
            ! write progress
            nloop = nloop + 1
            if (modulo(nloop,1000) == 0) write (error_unit,'(a)',advance='no') '.'
        end do

        write(error_unit,'(/,a,i0,a)') 'got ',nz,' zones'
        write(error_unit,'(a,i0,a,f5.3)') 'rejected ',nrejected,' zones: dP/P < ',delta_lgP
        close(unitno)

        call realloc_abuntime_arrays(nz,ierr)
        T = T*1.0e9
        
    contains
        subroutine realloc_abuntime_arrays(newsize,ierr)
            integer, intent(in) :: newsize
            integer, intent(out) :: ierr
            integer :: currentsize, oldsize
            real(dp), dimension(:), allocatable :: tmp_lgP
            real(dp), dimension(:,:), allocatable :: tmp_Y

            allocate(tmp_Y(nion,newsize),tmp_lgP(newsize), stat=ierr)
            if (ierr /= 0) return
            oldsize = size(lgP)

            currentsize = min(oldsize,newsize)
            tmp_lgP(1:currentsize) = lgP(1:currentsize)
            tmp_Y(:,1:currentsize) = Y(:,1:currentsize)
            deallocate(lgP,Y)
            allocate(lgP(newsize),Y(nion,newsize), stat=ierr)
            if (ierr /= 0) return
            lgP(1:newsize) = tmp_lgP(1:newsize)
            Y(:,1:newsize) = tmp_Y(:,1:newsize)

            deallocate(tmp_lgP,tmp_Y)
        end subroutine realloc_abuntime_arrays
    end subroutine read_abuntime

    subroutine reduce_abuntime(nz,nion,isos,lgP,Y,nnet,network,Yout,ierr)
        integer, intent(in) :: nz,nion
        character(len=iso_name_length), dimension(:), intent(in) :: isos
        real(dp), dimension(:), intent(in) :: lgP
        real(dp), dimension(:,:), intent(in) :: Y
        integer, intent(out) :: nnet
        character(len=iso_name_length), dimension(:), allocatable, intent(out) :: network
        real(dp), dimension(:,:), allocatable, intent(out) :: Yout
        integer, intent(out) :: ierr
        integer :: i
        logical, dimension(nion,nz) :: above_thresh
        logical, dimension(nion) :: ion_mask

        ierr = 0
        ion_mask = .TRUE.
        where (adjustl(isos) == 'n') ion_mask = .FALSE.
        do i = 1, nz
            above_thresh(:,i) = Y(:,i) > abundance_threshold*maxval(Y(:,i),ion_mask)
        end do
        
        ion_mask = .FALSE.        
        do i = 1, nion
            if (any(above_thresh(i,:))) ion_mask(i) = .TRUE.
        end do
        nnet = count(ion_mask)
        if (allocated(network)) deallocate(network)
        if (allocated(Yout)) deallocate(Yout)
        allocate(network(nnet),Yout(nnet,nz))
        network = pack(isos,ion_mask)
        do i = 1, nz
            Yout(:,i) = pack(Y(:,i),ion_mask)
        end do
    end subroutine reduce_abuntime

    subroutine expand_abuntime(nz,nion,isos,lgP,Y,lgP_increment, &
    & nzout,nnet,network,lgPout,Yout,ierr)
        use dStar_crust_def
        use hz90
        integer, intent(in) :: nz,nion
        character(len=iso_name_length), dimension(:), intent(in) :: isos
        real(dp), dimension(:), intent(in) :: lgP
        real(dp), dimension(:,:), intent(in) :: Y
        real(dp), intent(in) :: lgP_increment
        integer, intent(out) :: nzout, nnet
        character(len=iso_name_length), dimension(:), allocatable, intent(out) :: network
        real(dp), dimension(:), allocatable, intent(out) :: lgPout
        real(dp), dimension(:,:), allocatable, intent(out) :: Yout
        integer, intent(out) :: ierr
        character(len=iso_name_length), dimension(nion+HZ90_number) :: tmp_net
        integer :: indx, n_top, n_bottom, i
        real(dp) :: abuntime_lgP_min, abuntime_lgP_max
        real(dp), dimension(:,:), allocatable :: Y_HZ90
        
        ! make the network the union of the network and that of HZ90
        nnet = 0
        do i=1,nion
            nnet = nnet+1
            call integer_dict_define(isodir, adjustl(isos(i)), nnet, ierr)
            tmp_net(nnet) = isos(i)
        end do
        do i=1,HZ90_number
            call integer_dict_lookup(isodir, HZ90_network(i), indx, ierr)
            if (ierr == -1) then
                nnet = nnet+1
                call integer_dict_define(isodir, HZ90_network(i), nnet, ierr)
                tmp_net(nnet) = HZ90_network(i)
            end if
        end do
        allocate(network(nnet))
        network = tmp_net(1:nnet)
        
        abuntime_lgP_min = minval(lgP)
        abuntime_lgP_max = maxval(lgP)
        
        ! pad the table using HZ90
        n_top = (abuntime_lgP_min - crust_default_lgPmin)/lgP_increment
        n_bottom = (crust_default_lgPmax - abuntime_lgP_max)/lgP_increment
        nzout = n_top + nz + n_bottom
        allocate(lgPout(nzout),Yout(nnet,nzout))
        lgPout(1:n_top) = [ (crust_default_lgPmin+real(i-1)*lgP_increment, &
        & i=1,n_top) ]
        lgPout(n_top+1:n_top+nz) = lgP
        lgPout(n_top+nz+1:nzout) = [ (abuntime_lgP_max+real(i)*lgP_increment,i=1,n_bottom) ]

        Yout = 0.0_dp
        allocate(Y_HZ90(HZ90_number,n_top))
        call set_HZ90_composition(lgPout(1:n_top),Y_HZ90)
        do i = 1, HZ90_number
            call integer_dict_lookup(isodir,HZ90_network(i),indx,ierr)
            Yout(indx,1:n_top) = Y_HZ90(i,:)
        end do
!         Yout(1:nion,1:n_top) = spread(Y(1:nion,1),2,n_top)

        Yout(1:nion,n_top+1:n_top+nz) = Y(:,:)

        deallocate(Y_HZ90)
        allocate(Y_HZ90(HZ90_number,n_bottom))
        call set_HZ90_composition(lgPout(n_top+nz+1:nzout),Y_HZ90)
        do i = 1, HZ90_number
            call integer_dict_lookup(isodir,HZ90_network(i),indx,ierr)
            Yout(indx,n_top+nz+1:nzout) = Y_HZ90(i,:)
        end do
        
        ! smooth the transitions
        do i = 1, nnet
            where(lgPout <= abuntime_lgP_min .and.  &
                & lgPout >= abuntime_lgP_min-2*transition_width)
                Yout(i,:) =  Yout(i,n_top+1)* &
                &   cos(0.5*pi*(abuntime_lgP_min-lgPout)/2/transition_width) + &
                &   Yout(i,:)* (1.0- &
                &   cos(0.5*pi*(abuntime_lgP_min-lgPout)/2/transition_width))
            end where
            where(lgPout >= abuntime_lgP_max .and.  &
                & lgPout <= abuntime_lgP_max+2*transition_width)
                Yout(i,:) =  Yout(i,n_top+nz)* &
                &   cos(0.5*pi*(lgPout-abuntime_lgP_max)/2/transition_width) + &
                &   Yout(i,:)* (1.0- &
                &   cos(0.5*pi*(lgPout-abuntime_lgP_max)/2/transition_width))
            end where
        end do
    end subroutine expand_abuntime

    subroutine write_abuntime_cache(cache_filename,nz,nion,isos,T,lgP,Y,ierr)
        use nucchem_def, only: iso_name_length
        character(len=*), intent(in) :: cache_filename
        integer, intent(in) :: nz, nion
        character(len=iso_name_length), intent(in), dimension(:) :: isos ! nion
        real(dp), intent(in) :: T
        real(dp), intent(in), dimension(:) :: lgP   ! nz
        real(dp), intent(in), dimension(:,:) :: Y   ! nion, nz
        integer, intent(out) :: ierr
        integer :: unitno
    
        ierr = 0
        open(newunit=unitno, file=trim(cache_filename),action='write', &
        &   form='unformatted',iostat=ierr)
        if (failure('opening '//trim(cache_filename),ierr)) return

        write(unitno) nz
        write(unitno) nion
        write(unitno) isos
        write(unitno) T
        write(unitno) lgP
        write(unitno) Y
        close(unitno)
    end subroutine write_abuntime_cache

    subroutine read_abuntime_cache(cache_filename,nz,nion,isos, &
        &   T,lgP,Y,ierr)
        use nucchem_def, only: iso_name_length, composition_info_type

        character(len=*), intent(in) :: cache_filename
        integer, intent(out) :: nz,nion
        character(len=iso_name_length), intent(out), dimension(:), allocatable :: isos
        real(dp), intent(out) :: T
        real(dp), dimension(:), intent(out), allocatable :: lgP
        real(dp), dimension(:,:), intent(out), allocatable :: Y
        integer, intent(out) :: ierr
        integer :: unitno

        open(newunit=unitno,file=trim(cache_filename), &
        &   action='read',status='old',form='unformatted', iostat=ierr)
        if (failure('opening'//trim(cache_filename),ierr)) return

        read(unitno) nz
        read(unitno) nion
        allocate(isos(nion), lgP(nz), Y(nion,nz),stat=ierr)
        if (failure('allocating abuntime tables',ierr)) then
            close(unitno)
            return
        end if
        read(unitno) isos
        read(unitno) T
        read(unitno) lgP
        read(unitno) Y
        close(unitno)

    end subroutine read_abuntime_cache
      
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
