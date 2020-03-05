program test_P3_tab
    use iso_fortran_env, only: output_unit, error_unit
    use constants_lib
    use PPP_electron
    use num_lib, only: binary_search
    
    implicit none
    character(len=*), parameter :: datadir = '../../data/conductivity'
    type(electron_conductivity_tbl), pointer :: tab
    real(dp) :: Z, lgrho
    real(dp), dimension(5) :: lgT
    integer :: ierr, i, j
    
    lgT(:) = [ (real(j-1,dp)*0.5_dp + 7.0_dp, j=1,5) ]

    tab => PPP_tbl
    call load_PPP_electron_table(datadir,ierr)
    
    if (ierr /= 0) stop
    
    call construct_interpolation_coefficients(ierr)
    if (ierr /= 0) stop
    
    Z = 26.0_dp
    write(output_unit,'(tr11,5f8.3)') lgT(:)
    write(output_unit,'(tr11,40("-"))')
    do i = 1, 16
        lgrho = real(i-1,dp)/3.0_dp + 4.0_dp
        call do_one(10.0_dp**lgrho)
    end do
    
contains
    subroutine do_one(rho)
        real(dp), intent(in) :: rho
        real(dp) :: T, K
        real(dp), dimension(5) :: lgK
        integer :: j, ierr
        
        do j = 1, 5
            T = 10.0_dp**lgT(j)
            call eval_PPP_electron_table(rho,T,Z,K,ierr)
            if (ierr/= 0) stop
            lgK(j) = log10(K)
        end do
        write(output_unit,'(f8.3," | ",5f8.3)') log10(rho), lgK(:)
    end subroutine do_one
end program test_P3_tab
