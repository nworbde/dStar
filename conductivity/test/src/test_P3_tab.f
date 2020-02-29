program test_P3_tab
    use iso_fortran_env, only: output_unit, error_unit
    use constants_lib
    use num_lib, only: binary_search
    
    implicit none
    character(len=*), parameter :: datafile = '../data/condall06.dat'
    character(len=*), parameter :: minmaxfmt = '(a,f6.3,"; ",a,f6.3)'
    integer, parameter :: nZ = 15, nRho = 64, nT = 19
    real(dp), dimension(nZ) :: Zs
    real(dp), dimension(nRho) :: Rhos
    real(dp), dimension(nT) :: Ts
    real(dp), dimension(nT,nRho,nZ) :: lgK
    integer :: iounit, ierr
    integer :: iZ, iRho, lZ, lRho, lT
    real(dp) :: Ztest
    real(dp),dimension(1) :: Rhotest, Ttest, Ktest
    
    open(newunit=iounit,file=datafile,status='old',action='read',iostat=ierr)
    if (failure(ierr,'unable to open '//datafile)) stop
    
    write(error_unit,*) 'reading '//datafile//'...'
    read(iounit,*)  ! skip first line
    
    do iZ = 1, nZ
        read(iounit,*,iostat=ierr) Zs(iZ),Ts(:)
        if (failure(ierr,'unable to read Zs, Ts')) then
            write (error_unit,'(tr4,a,i2)') 'iZ = ',iZ
            stop
        end if
        do iRho = 1, nRho
            read(iounit,*,iostat=ierr) Rhos(iRho), lgK(:,iRho,iZ)
            if (failure(ierr,'unable to read Rho, lgK')) then
                write(error_unit,'(tr4,a,i2,", ",i2)') 'iZ, iRho = ',iZ,iRho
                stop
            end if
        end do
    end do
    close(iounit)
    write(error_unit,*) 'done'
    Zs = log10(Zs)
  
    write(output_unit,minmaxfmt) 'min(T) = ',minval(Ts),'max(T) = ',maxval(Ts)
    write(output_unit,minmaxfmt) 'min(Rho) = ',minval(Rhos), &
    &   'max(Rho) = ',maxval(Rhos)
    write(output_unit,minmaxfmt) 'min(K) = ',minval(lgK),'max(K) = ',maxval(lgK)
    
    Ztest = log10(26.5_dp)
    Rhotest = 6.1_dp
    Ttest = 7.5_dp
    
    lZ = binary_search(nZ,Zs,-1,Ztest)
    lRho = binary_search(nRho,Rhos,-1,Rhotest(1))
    lT = binary_search(nT,Ts,-1,Ttest(1))

    call interpolate(Ztest,Rhotest,Ttest,Ktest,ierr)
    
    write(output_unit,'(f7.3,tr7,f7.3)') lgK(lT:lT+1,lRho,lZ)
    write(output_unit,'(f7.3,tr7,f7.3)') lgK(lT:lT+1,lRho+1,lZ)

    write(output_unit,'(tr7,f7.3)') Ktest(1)

    write(output_unit,'(f7.3,tr7,f7.3)') lgK(lT:lT+1,lRho,lZ+1)
    write(output_unit,'(f7.3,tr7,f7.3)') lgK(lT:lT+1,lRho+1,lZ+1)    
    
contains
    subroutine interpolate(Z,Rho,T,K,ierr)
        use interp_2d_lib_db, only: interp_RGBI3P_db
        real(dp), intent(in) :: Z
        real(dp), intent(in), dimension(:) :: Rho, T
        real(dp), dimension(:), intent(out) :: K
        integer, intent(out) :: ierr
        real(dp), dimension(3,nT,nRho) :: wk
        real(dp), dimension(nT,nRho) :: lgK_Z
        integer, save :: loc = -1
        integer, parameter :: md = 1
        real(dp) :: wp, wm, norm
        integer :: dout
        
        loc = binary_search(nZ,Zs,loc,Z)
        wp = Z-Zs(loc); wm = Zs(loc+1)-Z
        norm = Zs(loc+1)-Zs(loc)
        lgK_Z = (wm*lgK(:,:,loc) + wp*lgK(:,:,loc+1))/norm
        
        dout = size(Rho)
        
        call interp_RGBI3P_db(md,nT,nRho,Ts,Rhos,lgK_Z,dout,T,Rho,K,ierr,wk)
    end subroutine interpolate
    
    function failure(ierr,msg)
        integer, intent(in) :: ierr
        character(len=*), intent(in) :: msg
        logical :: failure
        
        failure = .FALSE.
        if (ierr /= 0) then
            write(error_unit,*) 'FAILURE: '//msg
            failure = .TRUE.
        end if
    end function failure
end program test_P3_tab
