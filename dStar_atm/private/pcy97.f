! implements atmosphere model of Potekhin, Chabrier, and Yakovlev (1997).

module pcy97

	
	contains
	
	subroutine do_get_Teff(grav, Plight, Tb, Teff, flux)
		use constants_def
		real(dp), intent(in) :: grav	! surface gravity, in the local frame
		real(dp), intent(in) :: Plight	! pressure at which layer of light elements terminates
		real(dp), intent(in), dimension(:) :: Tb	! temperature at a base column
		real(dp), intent(out), dimension(:) :: Teff, flux	! effective temperature and flux
		real(dp) :: eta, g14
		real(dp), dimension(size(Tb)) :: Tb9, Teff6_4
		
		eta = Plight/1.193e34_dp
		g14 = grav * 1.0e-14_dp
		Tb9 = Tb*1.0e-9_dp
		call PCYfit(g14,eta,Tb9,Teff6_4)
		Teff = 1.0e6_dp*Teff6_4**0.25_dp
		flux = sigma_SB*Teff6_4*1.0e24_dp

		contains
		subroutine PCYfit(g14,eta,Tb9,Teff6_4)
			real(dp), intent(in) :: g14,eta
			real(dp), dimension(:), intent(in) :: Tb9
			real(dp), dimension(size(Tb9)), intent(out) :: Teff6_4	! (Teff/10^6 K)^4
			real(dp), dimension(size(Tb9)) :: a,Tstar,zeta, Te4Fe, Te4a
		
			Tstar = sqrt(7.0_dp*Tb9*sqrt(g14))
			zeta = Tb9-1.0e-3_dp*Tstar
			Te4Fe = g14*((7.0*zeta)**2.25_dp + (onethird*zeta)**1.25_dp)

			if (eta > 0) then
				Te4a = g14*(18.1_dp*Tb9)**2.42_dp
				a = tb9**fivethird * (1.2_dp + (5.3e-6_dp/eta)**0.38_dp)
				Teff6_4 = (a*Te4Fe + Te4a)/(a + 1.0_dp)
			else
				Teff6_4 = Te4Fe
			end if
		end subroutine PCYfit

	end subroutine do_get_Teff

end module pcy97
