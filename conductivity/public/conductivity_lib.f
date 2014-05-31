module conductivity_lib
	use conductivity_def

	contains
	
	subroutine get_thermal_conductivity(rho,T,chi,Gamma,eta,ionic,K,use_pcy,which_components)
		use nucchem_def, only: composition_info_type
		use eval_conductivity
		real(dp), intent(in) :: rho,T,Gamma,eta, chi
		type(composition_info_type), intent(in) :: ionic
		type(conductivity_components), intent(out) :: K
		logical, intent(in), optional :: use_pcy, which_components(num_conductivity_channels)
		integer :: which_ee
		logical, dimension(num_conductivity_channels) :: K_components

		which_ee = icond_sy06
		K_components = cond_use_all
		
		if (present(use_pcy)) then
			if (use_pcy) which_ee = icond_pcy
		end if
		
		if (present(which_components)) K_components = which_components
		
		call conductivity(rho,T,chi,Gamma,eta,ionic,K,which_ee,K_components)
	end subroutine get_thermal_conductivity

end module conductivity_lib
