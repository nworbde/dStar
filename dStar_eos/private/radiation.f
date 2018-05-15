module radiation_eos
    ! Implements eos for radiation
    
contains
    subroutine get_radiation_eos(rho,T,f,u,p,s,cv,dpr,dpt)
        use constants_def, only : dp, arad, onethird
        real(dp), intent(in) :: rho     ! g cm**-3
        real(dp), intent(in) :: T   ! K
        real(dp), intent(out) :: f, u, s    ! ergs/g
        real(dp), intent(out) :: p  ! dyn cm**-2
        real(dp), intent(out) :: cv ! ergs/K/g
        real(dp), intent(out) :: dpr,dpt    ! dp/dlnrho, dp/dlnT
        real(dp), parameter :: fourthird = 4.0_dp*onethird
        real(dp) :: aT3, aT4
        
        aT3 = arad*T**3
        aT4 = arad*T**4
        u = aT4 / rho
        p = onethird*aT4
        f = -onethird*aT4/rho
        s = fourthird*aT3/rho
        dpr = 0.0
        dpt = fourthird*aT4
        cv = 4.0*aT3/rho
    end subroutine get_radiation_eos

end module radiation_eos
