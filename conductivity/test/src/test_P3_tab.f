program test_P3_tab
    use iso_fortran_env, only: output_unit, error_unit
    use constants_lib
    use PPP_electron
    use num_lib, only: binary_search
    
    implicit none
    character(len=*), parameter :: datadir = '../../data/conductivity'
    character(len=*), parameter :: minmaxfmt = '(a,es10.3,"; ",a,es10.3)'
    type(electron_conductivity_tbl), pointer :: tab
    real(dp) :: rho, T, Z, K
    integer :: ierr
    
    tab => PPP_tbl
    call load_PPP_electron_table(datadir,ierr)
    
    if (ierr /= 0) stop
    
    write(output_unit,minmaxfmt)  &
    &   'min(T) = ',tab% Tmin,'max(T) = ',tab% Tmax
    write(output_unit,minmaxfmt)  &
    &   'min(rho) = ',tab% rhomin,'max(rho) = ',tab% rhomax
    write(output_unit,minmaxfmt)  &
    &   'min(K) = ',minval(tab% lgK),'max(K) = ',maxval(tab% lgK)
    
    call construct_interpolation_coefficients(ierr)
    if (ierr /= 0) stop
    
    rho = 2.0e6_dp
    T = 3.0e7_dp
    Z = 26.5_dp
    call eval_PPP_electron_table(rho,T,Z,K,ierr)
    
    write(output_unit,'(f4.1,tr1,3(f5.2,tr1))')  &
    &   Z, log10(rho), log10(T), log10(K)
end program test_P3_tab
