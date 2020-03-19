program test_P3_tab
    use iso_fortran_env, only: output_unit, error_unit
    use utils_lib, only: StrUpCase
    use exceptions_lib
    use constants_lib
    use nucchem_def
    use nucchem_lib
	use superfluid_def, only: max_number_sf_types, neutron_1S0
	use superfluid_lib
    use dStar_eos_lib
    use conductivity_lib
    use PPP_electron
    use num_lib, only: binary_search
    
    implicit none
    
    integer :: eos_handle,cond_handle,ierr,i,j,k
    integer :: Z(2), N(2)
    real(dp) :: A(2), X(2)
    integer, dimension(2) :: chem_ids, charged_ids
    real(dp), dimension(2) :: Y, Yion
    integer :: ncharged
    type(composition_info_type) :: ionic
    real(dp) :: chi, Xsum
    real(dp), dimension(num_dStar_eos_results) :: res
    character(len=*), parameter :: datadir = '../../data/conductivity'
    type(electron_conductivity_tbl), pointer :: tab
    real(dp) :: lgrho, lgT, rho, T, eta, Gamma, mu_e, K_e
    type(conductivity_components) :: kappa
    type(crust_eos_component), dimension(num_crust_eos_components) :: &
    &   eos_components
    integer :: phase
	real(dp), dimension(max_number_sf_types) :: Tcs
    character(len=iso_name_length) :: name(2)
    type(assertion) :: check_okay=assertion(scope='main')
    real(dp), dimension(13,5,11) :: diff
    integer :: loc(2), locz(3)
    
    call constants_init('',ierr)
    call check_okay% assert(ierr==0)
    call nucchem_init('../../data',ierr)
    call check_okay% assert(ierr==0)
    call dStar_eos_startup('../../data')
    call conductivity_startup('../../data')
    eos_handle = alloc_dStar_eos_handle(ierr)
    call check_okay% assert(ierr==0)
    cond_handle = alloc_conductivity_handle(ierr)
    call check_okay% assert(ierr==0)

    call conductivity_set_controls(cond_handle,which_ee_scattering=icond_sy06)
    
    Tcs = 1.0e9_dp
    Z = [ 2, 26 ]
    N = [ 2, 30 ]
    A = real(Z+N,dp)
    chem_ids = [ (get_nuclide_index_from_ZN(Z(k),N(k)),k=1,2) ]
    name = [ (nuclib% name(chem_ids(j)), j=1,2) ]
    do j = 1,2
        name(j)(1:1) = StrUpCase(name(j)(1:1))
    end do
    diff = 0.0_dp
    composition: do k = 1, 11
        X(1) = 0.1_dp*real(k-1,dp)
        X(2) = 1.0_dp -X(1)
        Y = X/A
        call compute_composition_moments(2,chem_ids,Y,ionic,Xsum, &
            & ncharged, charged_ids, Yion, exclude_neutrons=.TRUE.)
        
        write(output_unit,'(/,2("X(",a,") = ",f3.1,tr4))') &
        &   (trim(name(j)),X(j),j=1,2)
        write (output_unit, &
        &   '(3a6,6a11,a7/,3("======"),6("==========="),"=======")') &
            & 'lg(r)','lg(T)','<Z>','Gamma','eta_e', &
            & 'K_ee','K_ei','K_tot','K_table','diff'
            
        temperature: do j = 1, 5
            lgT = real(j-1,dp)*0.5_dp + 7.0_dp
            T = 10.0_dp**lgT
            
            density: do i = 1, 13
                lgrho = real(i-1,dp)/3.0_dp + 5.0_dp
                rho = 10.0_dp**lgrho
                chi = use_default_nuclear_size
                call eval_crust_eos(eos_handle,rho,T,ionic, &
                    & ncharged, charged_ids, Yion, Tcs,  &
                    & res, phase, chi, eos_components)
                eta = res(i_Theta) !1.0/TpT
                Gamma = res(i_Gamma)
                mu_e = res(i_mu_e)
                call get_thermal_conductivity(cond_handle,rho,T, &
                &   chi,Gamma,eta,mu_e,ionic,Tcs(neutron_1S0),kappa)
                call eval_PPP_electron_table(rho,T,sqrt(ionic% Z2),K_e,ierr)
                diff(i,j,k) = (K_e - kappa% electron_total)/kappa% electron_total
                write (output_unit, '(3f6.2,6es11.3,f7.3)') &
                    & lgrho,lgT,ionic%Z, &
                    & Gamma,mu_e*mev_to_ergs/boltzmann/T, &
                    & kappa%ee,kappa%ei,kappa%electron_total, K_e, diff(i,j,k)
            end do density
        end do temperature
        
        write(output_unit,'(a,f7.3)') 'rms(|diff|) = ',norm2(diff(:,:,k))/sqrt(real(13*5,dp))
        loc = maxloc(abs(diff(:,:,k)))
        write(output_unit,'(a,f7.3,2(a,f6.2))') 'max(|diff|) = ',maxval(abs(diff(:,:,k))), &
        &   ' at lg(rho) = ',real(loc(1)-1,dp)/3.0_dp + 5.0_dp, &
        &   '; lg(T) = ',real(loc(2)-1,dp)*0.5_dp + 7.0_dp
        loc = minloc(abs(diff(:,:,k)))
        write(output_unit,'(a,f7.3,2(a,f6.2))') 'min(|diff|) = ',minval(abs(diff(:,:,k))), &
        &   ' at lg(rho) = ',real(loc(1)-1,dp)/3.0_dp + 5.0_dp, &
        &   '; lg(T) = ',real(loc(2)-1,dp)*0.5_dp + 7.0_dp
        
    end do composition

    write(output_unit,'(/,a)') 'overall differences'
    write(output_unit,'(a,f7.3)') 'rms(|diff|) = ',norm2(diff)/sqrt(real(13*5*11,dp))
    locz = maxloc(abs(diff))
    write(output_unit,'(a,f7.3,3(a,f6.2))') 'max(|diff|) = ',maxval(abs(diff)), &
    &   ' at lg(rho) = ',real(locz(1)-1,dp)/3.0_dp + 5.0_dp, &
    &   '; lg(T) = ',real(locz(2)-1,dp)*0.5_dp + 7.0_dp, &
    &   '; X(He) = ',0.1_dp*real(locz(3)-1,dp)
    locz = minloc(abs(diff))
    write(output_unit,'(a,f7.3,3(a,f6.2))') 'min(|diff|) = ',minval(abs(diff)), &
    &   ' at lg(rho) = ',real(locz(1)-1,dp)/3.0_dp + 5.0_dp, &
    &   '; lg(T) = ',real(locz(2)-1,dp)*0.5_dp + 7.0_dp, &
    &   '; X(He) = ',0.1_dp*real(locz(3)-1,dp)
    
    call clear_composition(ionic)
    call conductivity_shutdown
    call dStar_eos_shutdown
    call nucchem_shutdown
    
end program test_P3_tab
