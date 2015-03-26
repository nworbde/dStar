module radiation_eos
    ! Implements eos for radiation
    
contains
    subroutine get_radiation_eos(rho,T,f,u,p,s,cv,dpr,dpt)
        use constants_def, only : dp, arad
		real(dp), intent(in) :: rho		! g cm**-3
		real(dp), intent(in) :: T	! K
		real(dp), intent(out) :: f, u, s	! ergs/g
		real(dp), intent(out) :: p	! dyn cm**-2
		real(dp), intent(out) :: cv	! ergs/K/g
        real(dp), intent(out) :: dpr,dpt    ! dp/dlnrho, dp/dlnT
        real(dp) :: aT4
        
        aT4 = arad*T**4
        u = aT4 / rho
        p = onethird*aT4
        f = -p/rho
        s = -4.0*f/T
        dpr = 0.0
        dpt = 4.0*p
        cv = 4.0*u/T
    end subroutine get_radiation_eos

end module radiation_eos
