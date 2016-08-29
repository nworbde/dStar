module neutron_eos
	! implements fit for pure neutron gas from Mackie & Baym 1977.
	! Should have a more up-to-date fmla.
	!
	! n is the local density of neutrons = Yn*n_b/(1-chi), where chi is the 
	! volume fraction occupied by the nucleus
	!
	! The MB77 fit is for a zero-temperature gas. To get the specific heat and 
	! entrozpy,we use the values for an ideal fermion gas. Not great, but need
	! to do something...
	! 

	contains
	subroutine MB77(n,T,Tns,f,u,p,s,cv,dpr,dpT)
!         use crust_eos_def
		use constants_def
		use fermi
		
		real(dp), intent(in) :: n		! cm**-3
		real(dp), intent(in) :: T	! K
		real(dp), intent(in) :: Tns	! neutron 1S0 critical temp [K]
		real(dp), intent(out) :: f, u, s	! ergs/neutron
		real(dp), intent(out) :: p	! dyn cm**-2
		real(dp), intent(out) :: cv	! ergs/K/neutron
		real(dp) :: dpr, dpT, mu_n	! dp/dlnN, dp/dlnT, chem. pot.
		real(dp), parameter, dimension(4) :: c = [1.2974,15.0298,-15.2343,7.4663]
!         real(dp), parameter :: kF_p = 0.0
		real(dp) :: k,pF,tau,v,R,dpdk,dkdn
		real(dp) :: lambda3, zeta	! zeta = chem. pot(ideal Fermi gas)/kT
		
		f = 0.0; u = 0.0; p = 0.0; s = 0.0; cv = 0.0; dpr = 0.0; dpT = 0.0
		if (n == 0.0) return
		
		k = (0.5*threepisquare*n)**onethird / cm_to_fm
		dkdn = onethird*k/n
! 		call sf_get_results(kF_p,k,Tc)
		
		u = k*(c(1) + k*(c(2) + k*(c(3) + k*c(4))))	! MeV per baryon
		u = u*mev_to_ergs
		p = onethird * n*k*(c(1) + k*(2.0*c(2)+k*(3.0*c(3)+k*4.0*c(4)))) ! Mev/cm**3
		p = p*mev_to_ergs
		dpdk = onethird*n*k*(c(1)+k*(4.0*c(2)+k*(9.0*c(3)+16.0*c(4)*k)))*mev_to_ergs
		dpr = (p/n + dpdk*dkdn)*n
		dpT = 0.0
		
		lambda3 = pi**2 * (hbar**2/(Mneutron*boltzmann*T))**1.5 / sqrt(2.0)
		zeta = ifermi12(n * lambda3)
		cv = 2.5*zfermi32(zeta)/zfermi12(zeta) - 4.5*zfermi12(zeta)/zfermim12(zeta)
		cv = cv * boltzmann
		
		s = fivethird*zfermi32(zeta)/zfermi12(zeta) - zeta
		s = s*boltzmann

		! Now set the superfluid reduction factor
! 		Tns = Tc(neutron_1S0)
		if (T < Tns) then
			tau = T/Tns
			v = sqrt(1.0-tau)*(1.456-0.157/sqrt(tau) + 1.764/tau)
			R = (0.4186+sqrt(1.007**2+(0.5010*v)**2))**2.5 * exp(1.456-sqrt(1.456**2+v**2))
		else
			R = 1.0
		end if

		cv = cv * R
		s = s*R		! ??? makes sense when very degenerate....
		
		f = u - T*s
	end subroutine MB77

end module neutron_eos
