module abuntime
    ! reader for abuntime files
    use const_def, only: dp, pi
    use utils_def
    use utils_lib
    use nucchem_def
    
    real(dp), parameter :: gravity = 1.85052e14_dp
    real(dp), parameter :: mdot = 8.8e4_dp*0.3_dp
    real(dp), parameter :: default_lgP_increment = 0.005_dp
    real(dp), parameter :: default_abundance_threshold = 0.01_dp
    real(dp), parameter :: default_lgPmin = 22.0_dp
    real(dp), parameter :: default_lgPmax = 33.5_dp
    
contains

    subroutine read_abuntime(filename, nz, nion, T, lgP, isos, Y, ierr, lgP_increment)
        use iso_fortran_env, only: error_unit, iostat_end
        use math_lib
        use exceptions_lib
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
        character(len=128) :: alert_msg
        type(alert) :: status = alert(scope='read_abuntime')
        type(failure) :: open_abuntime_failure = failure(scope='read_abuntime')
        type(failure) :: read_abuntime_failure = failure(scope='read_abuntime')
        type(failure) :: allocate_abuntime_failure = failure(scope='read_abuntime', &
            & message='allocating P, rho, Yion, ion_info, Xneut')

        delta_lgP = default_lgP_increment
        if (present(lgP_increment)) delta_lgP = lgP_increment

        open(newunit=unitno,file=trim(filename), &
        &   action='read',status='old', iostat=ierr)
        if (open_abuntime_failure% raised(ierr,message='opening'//trim(filename))) return

        read(unitno,'(1x,i5)') nion
        write(alert_msg,'(a,i5)') 'nion = ',nion
        call status% report(alert_msg)
        allocate(isos(nion))

        allocate(lgP(default_chunk_size), Y(nion,default_chunk_size), &
        &   stat=ierr)
        if (allocate_abuntime_failure% raised(ierr)) return

        nz = 1
        nloop = 1
        last_lgP = -100.0_dp
        nrejected = 0
        do
            ! size check
            if (nz > size(lgP)) then
                call realloc_abuntime_arrays(2*size(lgP),ierr)
                if (allocate_abuntime_failure% raised(ierr)) return
            end if
            read(unitno,'(1x,1e18.10,2e11.3,2X,2E18.10)',iostat=ios) &
                & time,T,rho,EFe,EFn
            if (ios == iostat_end) then
                nz = nz-1
                exit
            else if (read_abuntime_failure% raised (ios)) then
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
        write(error_unit,*)
        write(alert_msg,'(i0,a,i0,a)') nz,' zones accepted; ',nrejected,' zones rejected'
        call status% report(alert_msg)
        close(unitno)

        call realloc_abuntime_arrays(nz,ierr)
        if (allocate_abuntime_failure% raised(ierr)) return
        T = T*1.0e9_dp
        
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

    subroutine reduce_abuntime(nz,nion,isos,lgP,Y,nnet,network,Yout,ierr,abundance_threshold)
        use exceptions_lib
        integer, intent(in) :: nz,nion
        character(len=iso_name_length), dimension(:), intent(in) :: isos
        real(dp), dimension(:), intent(in) :: lgP
        real(dp), dimension(:,:), intent(in) :: Y
        integer, intent(out) :: nnet
        character(len=iso_name_length), dimension(:), allocatable, intent(out) :: network
        real(dp), dimension(:,:), allocatable, intent(out) :: Yout
        integer, intent(out) :: ierr
        real(dp), intent(in), optional :: abundance_threshold
        real(dp) :: min_abundance
        integer :: i
        logical, dimension(nion,nz) :: above_thresh
        logical, dimension(nion) :: ion_mask
        character(len=128) :: alert_msg
        type(alert) :: status = alert(scope='reduce_abuntime')

        min_abundance = default_abundance_threshold
        if (present(abundance_threshold)) min_abundance = abundance_threshold
        write(alert_msg,'(a,es9.2)')'setting min abundance to ',min_abundance
        call status% report(alert_msg)

        ierr = 0
        ion_mask = .TRUE.
        where (adjustl(isos) == 'n') ion_mask = .FALSE.
        do i = 1, nz
            above_thresh(:,i) = Y(:,i) > min_abundance*maxval(Y(:,i),ion_mask)
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

    subroutine extend_abuntime(nz,nion,isos,lgP,Y,lgP_increment,nzout,nnet,network,lgPout,Yout,ierr)
        use HZ90_comp
        use exceptions_lib
        implicit none
        real(dp), parameter :: lgP_blend_width = 0.2_dp
        integer, intent(in) :: nz,nion
        character(len=iso_name_length), dimension(nion), intent(in) :: isos
        real(dp), dimension(nz), intent(in) :: lgP
        real(dp), dimension(nion,nz), intent(in) :: Y
        real(dp), intent(in) :: lgP_increment
        integer, intent(out) :: nzout,nnet
        character(len=iso_name_length), dimension(:), allocatable, intent(out) :: network
        real(dp), dimension(:), allocatable, intent(out) :: lgPout
        real(dp), dimension(:,:), allocatable, intent(out) :: Yout
        integer, intent(out) :: ierr
        
        character(len=iso_name_length), dimension(:), allocatable :: tmp_net
        type(integer_dict), pointer :: isodir=>null()
        integer, dimension(:), allocatable :: indcs, indcsHZ90
        real(dp) :: lgPlow, lgPhigh
        real(dp), dimension(:,:), allocatable :: Yab,YHZ90
        real(dp), dimension(:), allocatable :: alpha, beta
        integer :: nlow, nhigh, i, indx
        type(assertion) :: dictionary_definition_okay = assertion(scope='extend_abuntime', &
        &   message='error in adding to network dictionary')
        character(len=128) :: alert_msg
        type(alert) :: status = alert(scope='extend_abuntime')
        
        ! construct union of abuntime network and HZ90
        allocate(tmp_net(nion+HZ90_number_isos))
        allocate(indcs(nion),indcsHZ90(HZ90_number_isos))
        nnet = 0
        do i = 1, nion
            nnet = nnet + 1
            call integer_dict_define(isodir,adjustl(isos(i)), nnet, ierr)
            call dictionary_definition_okay% assert(ierr==0)
            tmp_net(nnet) = isos(i)
            indcs(i) = nnet
        end do
        do i = 1, HZ90_number_isos
            call integer_dict_lookup(isodir, adjustl(HZ90_network(i)), indx, ierr)
            if (ierr == -1) then
                nnet = nnet + 1
                call integer_dict_define(isodir,adjustl(HZ90_network(i)), nnet, ierr)
                call dictionary_definition_okay% assert(ierr==0)
                tmp_net(nnet) = HZ90_network(i)
                indcsHZ90(i) = nnet
            else
                indcsHZ90(i) = indx
            end if
        end do
        allocate(network(nnet))
        network = tmp_net(1:nnet)
        deallocate(tmp_net)
        
        ! extend lg P to lower, higher values
        lgPlow = minval(lgP)
        lgPhigh = maxval(lgP)
        nlow = max(int((lgPlow - default_lgPmin)/lgP_increment),1)
        nhigh = max(int((default_lgPmax - lgPhigh)/lgP_increment),1)
        nzout = nlow+nz+nhigh
        allocate(lgPout(nzout),alpha(nzout),beta(nzout))
        allocate(Yab(nion,nzout),YHZ90(HZ90_number_isos,nzout))

        lgPout(1:nlow) = (/(minval(lgP) - lgP_increment*i,i=nlow,1,-1)/)
        lgPout(nlow+1:nlow+nz) = lgP
        lgPout(nlow+nz+1:nzout) = (/(maxval(lgP) + lgP_increment*i, i= 1, nhigh)/)
        
        write(alert_msg,'(2(a,f7.3,"--",f7.3))') 'extending lgP range from ',lgPlow,lgPhigh, &
        &   ' to ',minval(lgPout),maxval(lgPout)
        call status% report(alert_msg)

        Yab(:,1:nlow) = spread(Y(:,1),2,nlow)
        Yab(:,nlow+1:nlow+nz) = Y
        Yab(:,nlow+nz+1:nzout) = spread(Y(:,nz),2,nhigh)        
        call do_generate_HZ90_table(lgPout,YHZ90)
        
        allocate(Yout(nnet,nzout),source = 0.0_dp)
        where(lgPout <= lgPhigh)
            beta = 0.0_dp
        else where(lgPout > lgPhigh+lgP_blend_width)
            beta = 1.0_dp
        else where(lgPout > lgPhigh .and. lgPout <= lgPhigh+lgP_blend_width)
            beta = (lgPout(:)-lgPhigh)/lgP_blend_width
        end where
        alpha = 1.0_dp-beta
        do i = 1, nion
            Yout(indcs(i),:) = Yab(i,:)*alpha(:)
        end do
        do i = 1, HZ90_number_isos
            Yout(indcsHZ90(i),:) = Yout(indcsHZ90(i),:) + YHZ90(i,:)*beta(:)
        end do
        
        deallocate(YHZ90,Yab,indcs,indcsHZ90,alpha,beta)
        
    end subroutine extend_abuntime
    
end module abuntime
