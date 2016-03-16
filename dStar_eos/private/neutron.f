module neutron_eos
	! implements fit for pure neutron gas.  There is the Mackie & Baym (1977) 
	! formalism as well as the Skyrme parameter sets from B.A. Brown.
	! 
	! A wrapper routine selects which one to call 
	!
	! n is the local density of neutrons = Yn*n_b/(1-chi), where chi is the 
	! volume fraction occupied by the nucleus
	!
	! The MB77 fit is for a zero-temperature gas. To get the specific heat and 
	! entropy, we use the values for an ideal fermion gas. Not great, but need 
	! to do something...
	! 

contains
		
	subroutine get_neutron_eos(rq,n,T,Tns,f,u,p,s,cv,dpr,dpt)
		use constants_def
		use fermi
		use superfluid_def, only: max_number_sf_types, neutron_1S0
		use superfluid_lib, only: sf_get_results
		use dStar_eos_private_def
		type(dStar_eos_general_info), pointer :: rq		
		real(dp), intent(in) :: n		! cm**-3
		real(dp), intent(in) :: T	! K
		real(dp), intent(in) :: Tns	! neutron 1S0 critical temp [K]
		real(dp), intent(out) :: f, u, s	! ergs/neutron
		real(dp), intent(out) :: p	! dyn cm**-2
		real(dp), intent(out) :: cv	! ergs/K/neutron
		real(dp), intent(out) :: dpr, dpT ! dp/dlnN, dp/dlnT
        real(dp), parameter :: kF_p = 0.0
		real(dp) :: k,tau,v,R,meff
		real(dp) :: lambda3, zeta	! zeta = chem. pot(ideal Fermi gas)/kT

		f = 0.0; u = 0.0; p = 0.0; s = 0.0; cv = 0.0; dpr = 0.0; dpT = 0.0
		if (n == 0.0) return
		
        ! get wavenumber in inverse fm
        k = (0.5*threepisquare*n)**onethird / cm_to_fm
		meff = 1.0_dp
		if (rq% use_skyrme_for_neutrons) then
			call get_skyrme(n,T,u,p,dpr,dpT,meff)
		else
			call get_MB77(n,T,u,p,dpr,dpT,wavenumber=k)
		end if
		
		! correct (approximately) for non-degeneracy
		lambda3 = pi**2 * (hbar**2/(Mneutron*meff*boltzmann*T))**1.5 / sqrt(2.0)
		zeta = ifermi12(n * lambda3)
		cv = 2.5*zfermi32(zeta)/zfermi12(zeta)  &
			&	- 4.5*zfermi12(zeta)/zfermim12(zeta)
		cv = cv * boltzmann
		
		s = fivethird*zfermi32(zeta)/zfermi12(zeta) - zeta
		s = s*boltzmann

		if (T < Tns) then
			tau = T/Tns
			v = sqrt(1.0-tau)*(1.456-0.157/sqrt(tau) + 1.764/tau)
			R = (0.4186+sqrt(1.007**2+(0.5010*v)**2))**2.5 * &
				 &	exp(1.456-sqrt(1.456**2+v**2))
		else
			R = 1.0
		end if

		cv = cv*R
		s = s*R
		
		f = u - T*s		
	end subroutine get_neutron_eos

	subroutine get_skyrme(n,T,u,p,dpr,dpt,meff)
		use constants_def
		use dStar_eos_def
		use skyrme
		real(dp), intent(in) :: n	! cm**-3
		real(dp), intent(in) :: T   ! K
		real(dp), intent(out) :: u	! ergs/neutron
		real(dp), intent(out) :: p  ! dyn cm**-2
		real(dp), intent(out) :: dpr,dpt,meff	! dp/dlnN, dp/dlnT, m*/m
		real(dp) :: n_n,dur
		
		u = 0.0; p = 0.0; dpr = 0.0; dpT = 0.0
		
		! convert to nuclear units
		n_n = n / density_n
		
		call eval_skyrme_eos(skyrme_neutron,n_n,u,p,dur,dpr,meff)
		
		! convert back to cgs
		u = u*mev_to_ergs
		p = p*pressure_n
		dur = dur*mev_to_ergs
		dpr = dpr*pressure_n		
	end subroutine get_skyrme

	subroutine get_MB77(n,T,u,p,dpr,dpT,wavenumber)
		use constants_def
		use fermi
		
		real(dp), intent(in) :: n		! cm**-3
		real(dp), intent(in),optional :: wavenumber		! fm**-1; will be calculated if not provided
		real(dp), intent(in) :: T	! K
		real(dp), intent(out) :: u	! ergs/neutron
		real(dp), intent(out) :: p	! dyn cm**-2
		real(dp), intent(out) :: dpr, dpT ! dp/dlnN, dp/dlnT
		real(dp) :: mu_n	! chem. pot.
		real(dp), parameter, dimension(4) :: c = [1.2974,15.0298,-15.2343,7.4663]
		real(dp) :: k,dpdk,dkdn
		
		u = 0.0; p = 0.0; dpr = 0.0; dpT = 0.0
		if (n == 0.0) return
		
		if (present(wavenumber)) then
			k = wavenumber
		else
			k = (0.5*threepisquare*n)**onethird / cm_to_fm
		end if

		dkdn = onethird*k/n		
		u = k*(c(1) + k*(c(2) + k*(c(3) + k*c(4))))	! MeV per baryon
		u = u*mev_to_ergs
		p = onethird * n*k*(c(1) + k*(2.0*c(2)+k*(3.0*c(3)+k*4.0*c(4)))) ! Mev/cm**3
		p = p*mev_to_ergs
		dpdk = onethird*n*k*(c(1)+k*(4.0*c(2)+k*(9.0*c(3)+16.0*c(4)*k)))*mev_to_ergs
		dpr = (p/n + dpdk*dkdn)*n
		dpT = 0.0
		
	end subroutine get_MB77

end module neutron_eos
