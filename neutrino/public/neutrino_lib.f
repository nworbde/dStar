module neutrino_lib
use neutrino_def
use superfluid_def, only: max_number_sf_types

contains

subroutine get_core_neutrino_emissivity(nn,np,T,Tcs,use_mode,epsilon_nu)
    use eval_core_neutrino
    real(dp), intent(in) :: nn,np,T,Tcs(max_number_sf_types)
    logical, intent(in), dimension(num_core_nu_channels) :: use_mode
    type(core_neutrino_emissivity_channels), intent(out) :: epsilon_nu

    call emissivity(nn,np,T,Tcs,use_mode,epsilon_nu)
end subroutine get_core_neutrino_emissivity

subroutine get_crust_neutrino_emissivity(rho,T,ionic,chi,Tcn,epsilon,use_mode)
    use nucchem_def, only: composition_info_type
    use eval_crust_neutrino
    real(dp), intent(in) :: rho,T,chi
    type(composition_info_type), intent(in) :: ionic
    real(dp), intent(in) :: Tcn ! neutron critical tempererature
    type(crust_neutrino_emissivity_channels), intent(out) :: epsilon
    logical, intent(in), dimension(num_crust_nu_channels), optional :: use_mode
    logical, dimension(num_crust_nu_channels) :: channels
    
    channels = crust_nu_full_cooling
    if (present(use_mode)) channels = use_mode
    
    call emissivity(rho,T,ionic,chi,Tcn,channels,epsilon)
end subroutine get_crust_neutrino_emissivity

end module neutrino_lib
