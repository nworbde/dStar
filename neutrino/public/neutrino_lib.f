module core_neutrino_lib
use core_neutrino_def
use superfluid_def, only: max_number_sf_types

contains

subroutine get_core_neutrino_emissivity(nn,np,T,Tcs,use_mode,epsilon_nu)
	use eval_core_neutrino
	real(dp), intent(in) :: nn,np,T,Tcs(max_number_sf_types)
	logical, intent(in), dimension(num_core_nu_channels) :: use_mode
	type(core_neutrino_emissivity_channels), intent(out) :: epsilon_nu

	call emissivity(nn,np,T,Tcs,use_mode,epsilon_nu)
end subroutine get_core_neutrino_emissivity

end module core_neutrino_lib
