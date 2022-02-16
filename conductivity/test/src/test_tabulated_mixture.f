program test_tabulated_mixture
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
    
    implicit none
    integer, parameter :: NX=11,Nrho=13,NT=5
    real(dp), parameter :: Xmin=0.0_dp, Xmax = 1.0_dp
    real(dp), parameter :: lgrho_min=5.0_dp, lgrho_max = 9.0_dp
    real(dp), parameter :: lgT_min=7.0_dp, lgT_max = 9.0_dp
    integer :: eos_handle,cond_handle,ierr,i,j,k
    integer :: Z(2), N(2)
    real(dp) :: A(2), X(2)
    integer, dimension(2) :: chem_ids, charged_ids
    real(dp), dimension(2) :: Y, Yion
    integer :: ncharged
    type(composition_info_type) :: ionic
    real(dp) :: chi, Xsum
    real(dp), dimension(num_dStar_eos_results) :: res
    character(len=*), parameter :: datadir = '../../data'
    real(dp), dimension(NX) :: X1
    real(dp), dimension(Nrho) :: lgrho
    real(dp), dimension(NT) :: lgT
    real(dp) :: rho,T,Gamma,eta,mu_e,K_e,K_e1,K_e2,Z1,Z2,A1,A2
    type(conductivity_components) :: kappa
    type(crust_eos_component), dimension(num_crust_eos_components) :: eos_components
    integer :: phase
	real(dp), dimension(max_number_sf_types) :: Tcs
    character(len=iso_name_length) :: name(2)
    type(assertion) :: check_okay=assertion(scope='main')
    real(dp), dimension(Nrho,NT,NX) :: diff
    integer :: loc(2), locz(3)
    character(len=16) :: method ! 'mean Z', 'rms Z', 'Z2/Z', 'density averaged'
    
    method='mean Z'
    if (command_argument_count() > 0) then
        call get_command_argument(1,method)
    end if
    
    ! turn off all status message
    call set_verbosity(1)
    
    ! initialize microphysics
    call constants_init('',ierr)
    call check_okay% assert(ierr==0)
    call nucchem_init(datadir,ierr)
    call check_okay% assert(ierr==0)
    call dStar_eos_startup(datadir)
    eos_handle = alloc_dStar_eos_handle(ierr)
    call check_okay% assert(ierr==0)
    ! Suppress the warnings for degenerate ions
    call dStar_eos_set_controls(eos_handle,suppress_warnings=.TRUE.)
    call conductivity_startup(datadir)
    cond_handle = alloc_conductivity_handle(ierr)
    call check_okay% assert(ierr==0)
    call conductivity_set_controls(cond_handle,which_ee_scattering=icond_pcy)    

    diff = 0.0_dp
    ! set superfluid temps to a uniform value; we aren't using these
    Tcs = 1.0e9_dp

    ! set composition
    Z = [ 26, 6 ]
    N = [ 30, 6 ]
    A = real(Z+N,dp)
    chem_ids = [ (get_nuclide_index_from_ZN(Z(k),N(k)),k=1,2) ]
    name = [ (nuclib% name(chem_ids(j)), j=1,2) ]
    name(1:2)(1:1) = [ (StrUpCase(name(j)(1:1)), j=1,2) ]
    Z1 = real(Z(1),dp); Z2 = real(Z(2),dp); A1 = real(A(1),dp); A2 = real(A(2),dp)

    ! set grid axes
    X1 = linspace(NX,Xmin,Xmax)
    lgrho = linspace(Nrho,lgrho_min,lgrho_max)
    lgT = linspace(NT,lgT_min,lgT_max)
    
    write(output_unit,'(a)') 'method: '//trim(method)
    composition: do k = 1, NX
        X(1) = X1(k)
        X(2) = 1.0_dp -X(1)
        Y = X/A
        call compute_composition_moments(2,chem_ids,Y,ionic,Xsum, &
            & ncharged, charged_ids, Yion, exclude_neutrons=.TRUE.)
        
        ! table header
        write(output_unit,'(2("X(",a,") = ",f3.1,tr4))') &
        &   (trim(name(j)),X(j),j=1,2)
        write(output_unit,'(2(a,f6.2,tr4))') '<Z> = ',ionic%Z,'<A> = ',ionic% A
        write (output_unit, &
        &   '(2a6,7a11,a7/,2("======"),7("==========="),"=======")') &
            & 'lg(r)','lg(T)','Gamma','eta_e', &
            & 'K_ee','K_ei','K_eQ','K_tot','K_table','diff'
            
        temperature: do j = 1, NT
            T = 10.0_dp**lgT(j)
            density: do i = 1, Nrho
                rho = 10.0_dp**lgrho(i)
                chi = use_default_nuclear_size
                call eval_crust_eos(eos_handle,rho,T,ionic, &
                    & ncharged, charged_ids, Yion, Tcs,  &
                    & res, phase, chi, eos_components)
                eta = res(i_Theta) !1.0/TpT
                Gamma = res(i_Gamma)
                mu_e = res(i_mu_e)
                call get_thermal_conductivity(cond_handle,rho,T, &
                &   chi,Gamma,eta,mu_e,ionic,Tcs(neutron_1S0),kappa)
                K_e = get_tabulated_conductivity(trim(method),rho,T,Y,Z1,Z2,A1,A2,ionic)
                if (kappa%eQ > 0.0_dp) K_e = 1.0/(1.0/K_e + 1.0/kappa%eQ)
                diff(i,j,k) = (K_e - kappa% electron_total)/kappa% electron_total
                write (output_unit, '(2f6.2,7es11.3,f7.3)') &
                    & lgrho(i),lgT(j), &
                    & Gamma,mu_e*mev_to_ergs/boltzmann/T, &
                    & kappa%ee,kappa%ei,kappa%eQ,kappa%electron_total, &
                    & K_e, diff(i,j,k)
            end do density
        end do temperature
        
        ! comparison for this composition
        write(output_unit,'(a,f7.3)') 'rms(|diff|) = ',norm2(diff(:,:,k))/sqrt(real(Nrho*5,dp))
        loc = maxloc(abs(diff(:,:,k)))
        write(output_unit,'(a,f7.3,2(a,f6.2))') 'max(|diff|) = ',maxval(abs(diff(:,:,k))), &
        &   ' at lg(rho) = ',lgrho(loc(1)), &
        &   '; lg(T) = ',lgT(loc(2))
        loc = minloc(abs(diff(:,:,k)))
        write(output_unit,'(a,f7.3,2(a,f6.2))') 'min(|diff|) = ',minval(abs(diff(:,:,k))), &
        &   ' at lg(rho) = ',lgrho(loc(1)), &
        &   '; lg(T) = ',lgT(loc(2))
        write(output_unit,*)
    end do composition

    ! global comparison
    write(output_unit,'(a)') 'overall differences'
    write(output_unit,'(a,f7.3)') 'rms(|diff|) = ',norm2(diff)/sqrt(real(Nrho*NT*NX,dp))
    locz = maxloc(abs(diff))
    write(output_unit,'(a,f7.3,2(a,f6.2),a,a,a,f6.2)') 'max(|diff|) = ',maxval(abs(diff)), &
    &   ' at lg(rho) = ',lgrho(locz(1)), &
    &   '; lg(T) = ',lgT(locz(2)), &
    &   '; X(',trim(name(1)),') = ',X1(locz(3))
    locz = minloc(abs(diff))
    write(output_unit,'(a,f7.3,2(a,f6.2),a,a,a,f6.2)') 'min(|diff|) = ',minval(abs(diff)), &
    &   ' at lg(rho) = ',lgrho(locz(1)), &
    &   '; lg(T) = ',lgT(locz(2)), &
    &   '; X(',trim(name(1)),') = ',X1(locz(3))

    call conductivity_shutdown
    call dStar_eos_shutdown
    call nucchem_shutdown
    
contains
    function get_tabulated_conductivity(method,rho,T,Y,Z1,Z2,A1,A2,ionic) result(K_e)
        character(len=*), intent(in) :: method
        real(dp), intent(in) :: rho,T,Y(2),Z1,Z2,A1,A2
        type(composition_info_type), intent(in) :: ionic
        integer :: ierr
        real(dp) :: K_e, K_e1, K_e2, rho1, rho2
        type(assertion) :: interp_okay=assertion(scope='get_tabulated_conductivity')

        select case(method)
        case ('mean Z')
            call eval_PPP_electron_table(rho,T,ionic% Z,K_e,ierr)
            call interp_okay% assert(ierr == 0)
        case ('rms Z')
            call eval_PPP_electron_table(rho,T,sqrt(ionic% Z2),K_e,ierr)
            call interp_okay% assert(ierr == 0)
        case ('Z2/Z')
            call eval_PPP_electron_table(rho,T,ionic% Z2/ionic% Z,K_e,ierr)
            call interp_okay% assert(ierr == 0)
        case ('density averaged')
            ! Average over the electron conductivity, assuming that ei scattering is 
            ! proportional to n_i = Y_i*rho/m_u = density of each nuclide. The density 
            ! is adjusted so that n_e is the same for each call.
            rho1 = ionic% Ye*A1/Z1 * rho
            call eval_PPP_electron_table(rho1,T,Z1,K_e1,ierr)
            call interp_okay% assert(ierr == 0)
            rho2 = ionic% Ye*A2/Z2 * rho
            call eval_PPP_electron_table(rho2,T,Z2,K_e2,ierr)
            call interp_okay% assert(ierr == 0)
            K_e = ionic% Ye/(Z1*Y(1)/K_e1 + Z2*Y(2)/K_e2)
        case default
            stop 'bad selection'
        end select
    end function get_tabulated_conductivity
    function linspace(jmax,xmin,xmax)
        integer, intent(in) :: jmax
        real(dp), intent(in) :: xmin, xmax
        real(dp), dimension(jmax) :: linspace
        linspace = [ (xmin + (xmax-xmin)*real(j-1,dp)/real(jmax-1,dp), j=1,jmax) ]
    end function linspace
end program test_tabulated_mixture
