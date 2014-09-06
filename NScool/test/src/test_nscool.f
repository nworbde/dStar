program test_nscool
    use constants_def, only: dp
    use constants_lib
    use nucchem_def
    use nucchem_lib
    use composition_models
    
    integer, parameter :: Ntab = 500
    integer :: i, ierr
    real(dp), dimension(Ntab) :: lgP
    real(dp), dimension(:,:), allocatable :: X
    
    call constants_init('',ierr)
    call nucchem_init('../../data',ierr)
    
    if (ierr /= 0) stop
    
    lgP = [(26.5+5.0*real(i-1,dp)/real(Ntab-1,dp),i=1,Ntab)]
    call HZ90(lgP,X)
    
    write(*,'(21a9,/)') 'lg P',HZ90_network(:),'sum'
    do i = 1,Ntab
        write(*,'(21f9.6)') lgP(i),X(:,i),sum(X(:,i))
        if (mod(i,20) == 0) then
            write(*,'(21a9,/)') 'lg P',HZ90_network(:),'sum'
        end if            
    end do
    
    call nucchem_shutdown

end program test_nscool
