program test_nscool
    use constants_def, only: dp
    use constants_lib
    use nucchem_def
    use nucchem_lib
    use composition_models
    
    integer, parameter :: Ntab = 500
    integer :: i, ierr
    real(dp), dimension(Ntab) :: lgrho
    real(dp), dimension(:,:), allocatable :: X
    
    call constants_init('',ierr)
    call nucchem_init('../../data',ierr)
    
    if (ierr /= 0) stop
    
    lgrho = [(8.5+5.0*real(i-1,dp)/real(Ntab-1,dp),i=1,Ntab)]
    call HZ90(lgrho,X)
    
    write(*,'(20a9,/)') 'lg rho',HZ90_network(:)
    do i = 1,Ntab
        write(*,'(20f9.6)') lgrho(i),X(:,i)
        if (mod(i,20) == 0) then
            write(*,'(20a9,/)') 'lg rho',HZ90_network(:)
        end if            
    end do
    
    call nucchem_shutdown

end program test_nscool
