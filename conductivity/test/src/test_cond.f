program test_cond
    use, intrinsic :: iso_fortran_env, only: output_unit
    use constants_lib
    use nucchem_def
    use nucchem_lib
    use dStar_eos_lib
    use conductivity_lib
    
    integer :: eos_handle,cond_handle,ierr,i,j
    ! composition taken from HZ090 for Fe-chain accreted crust
    integer, dimension(2:14), parameter ::  zz = [2,2,2,2,26,26,26,26,24,20,14,24,24], &
    &   aa = [4,4,4,4,56,56,56,56,56,56,46,96,96]
    real(dp), dimension(2:14), parameter :: xn = [0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.07,0.76,0.80]
    integer, dimension(2) :: Z,  N, chem_ids, charged_ids
    real(dp), dimension(2) :: Y, Yion
    integer :: ncharged
    type(composition_info_type) :: ionic
    real(dp) :: chi,Xsum
    real(dp), dimension(num_dStar_eos_results) :: res
    
    call constants_init('',ierr)
    call nucchem_init('../../data',ierr)
    call dStar_eos_startup('../../data')
    call conductivity_startup('../../data')
    eos_handle = alloc_dStar_eos_handle(ierr)
    cond_handle = alloc_conductivity_handle(ierr)
    call conductivity_set_controls(cond_handle, &
    &   include_electrons=.TRUE., &
    &   include_neutrons=.TRUE., &
    &   include_superfluid_phonons=.TRUE., &
    &   include_photons=.TRUE.)
    
    write (output_unit, '(5a6,10a14,/,5("======"),10("=============="))') &
        & 'lg(r)','lg(T)','<Z>','<A>','Y_n', &
        & 'Gamma','eta_e', &
        & 'K_tot','K_ee','K_ei','K_eQ', &
        & 'K_sF','K_nQ', 'K_n,phn', &
        & 'kappa_rad'
        
    do i = 2,14
        N = [1,aa(i)-zz(i)]
        Z = [0,zz(i)]
        chem_ids = [(get_nuclide_index_from_ZN(Z(j),N(j)),j=1,2)]

        Y = [xn(i),(1.0-xn(i))/real(aa(i),dp)]
        call compute_composition_moments(2,chem_ids,Y,ionic,Xsum, &
            & ncharged, charged_ids, Yion, exclude_neutrons=.TRUE.)
        ionic% Q = 4.0
        call do_one(10.0_dp**i)
    end do
    call clear_composition(ionic)
    call conductivity_shutdown()
    call dStar_eos_shutdown()
    call nucchem_shutdown

    contains
    subroutine do_one(rho)
		use superfluid_def, only: max_number_sf_types, neutron_1S0
		use superfluid_lib
        real(dp), intent(in) :: rho
        real(dp) :: Gamma,eta,f,u,p,s,cv,chi_rho,chi_T
        integer :: phase
        real(dp) :: K, lgr, lgT, T, TpT, mu_e
        type(conductivity_components) :: kappa
        type(crust_eos_component), dimension(num_crust_eos_components) :: eos_components
        integer :: ii
		real(dp), dimension(max_number_sf_types) :: Tcs
        
        ! for checking superfluid reduction/turn-on of phonon cond.
        Tcs = 1.0e8_dp

        lgr = log10(rho)
        chi = use_default_nuclear_size
        do ii = 1, 16
            lgT = 7.0_dp + real(ii-1,dp)/10.0_dp
            T = 10.0**lgT
            call eval_crust_eos(eos_handle,rho,T,ionic, &
                & ncharged, charged_ids, Yion, Tcs, res, phase, chi, eos_components)
            eta = res(i_Theta) !1.0/TpT
            Gamma = res(i_Gamma)
            mu_e = res(i_mu_e)
            call get_thermal_conductivity(cond_handle,rho,T, &
            &   chi,Gamma,eta,mu_e,ionic,Tcs(neutron_1S0),kappa)
            write (output_unit, '(5f6.2,10es14.6)') &
                & lgr,lgT,ionic%Z,ionic%A,ionic%Yn, &
                & Gamma,mu_e*mev_to_ergs/boltzmann/T, &
                & kappa% total,kappa%ee,kappa%ei,kappa%eQ, &
                & kappa% sf, kappa%nQ, kappa%np, kappa% kap
        end do
    end subroutine do_one
end program test_cond
