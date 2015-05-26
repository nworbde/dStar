module conductivity_lib
    use conductivity_def

contains

    subroutine get_thermal_conductivity( &
    &   rho,T,chi,Gamma,eta,mu_e,ionic,K,use_pcy,use_page,which_components)
        use nucchem_def, only: composition_info_type
        use eval_conductivity
        real(dp), intent(in) :: rho,T,Gamma,eta, mu_e, chi
        type(composition_info_type), intent(in) :: ionic
        type(conductivity_components), intent(out) :: K
        logical, intent(in), optional :: use_pcy, use_page, which_components(num_conductivity_channels)
        integer :: which_ee, which_eQ
        logical, dimension(num_conductivity_channels) :: K_components

        which_ee = icond_sy06
        which_eQ = icond_eQ_potekhin
        K_components = cond_use_all

        if (present(use_pcy)) then
            if (use_pcy) which_ee = icond_pcy
        end if

        if (present(use_page)) then
            if (use_page) which_eQ = icond_eQ_page
        end if

        if (present(which_components)) K_components = which_components

        call conductivity(rho,T,chi,Gamma,eta,mu_e,ionic,K,which_ee,which_eQ,K_components)
    end subroutine get_thermal_conductivity

    subroutine get_core_thermal_conductivity(nn,np,mneff,mpeff,T,Tcs,K)
        ! nn  := neutron density/fm**-3
        ! np  := proton density/fm**-3
        ! mneff := neutron effective mass/neutron rest mass
        ! mpeff := proton effective mass/proton rest mass
        ! T   := temperature (K)
        ! Tcs := critical temperatures (K)
        ! K   := thermal conductivity (cgs)
        use eval_core_conductivity
        use superfluid_def, only: max_number_sf_types
        real(dp), intent(in) :: nn,np,mneff,mpeff,T,Tcs(max_number_sf_types)
        real(dp), intent(out) :: K

        K = neutron_conductivity(nn,np,mneff,mpeff,T,Tcs)
    end subroutine get_core_thermal_conductivity

end module conductivity_lib
