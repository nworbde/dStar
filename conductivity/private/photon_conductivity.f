module photon_conductivity
    use conductivity_def
    implicit none
    
contains
    
    function electron_scattering(eta_e,theta,Ye) result(kTh)
        use constants_def
        real(dp), intent(in) :: eta_e   ! electron chemical pot./kT
        real(dp), intent(in) :: theta ! T/mc**2
        real(dp), intent(in) :: Ye
        real(dp) :: kTh
        real(dp) :: xi, xi2, theta2
        real(dp) :: t1, t2, t3, Gbar
    
        xi  = exp(0.8168*eta_e - 0.05522*eta_e**2)
        xi2 = xi**2
        theta2 = theta**2
        t1  = 1.129 + 0.2965*xi - 0.005594*xi2
        t2  = 11.47 + 0.3570*xi + 0.1078*xi2
        t3  = -3.249 + 0.1678*xi - 0.04706*xi2
  
        Gbar    =  t1 + t2*theta + t3*theta2
        kTh = Thomson*avogadro*Ye/Gbar
    end function electron_scattering

    ! Calculates the free-free Gaunt factor using a fitting
    ! formula described in Schatz et al. (1999)
    ! Agrees to 10% with Itoh et al. 1991
    !
    function freefree(rho,T,eta_e,ionic) result(kff)
        use nucchem_def, only: composition_info_type
        use constants_def
        real(dp), intent(in) :: rho,T,eta_e
        type(composition_info_type), intent(in) :: ionic
        real(dp) :: kff
    
        real(dp), parameter :: bfac = 0.753
        real(dp) :: rY,T8,T8_32
        real(dp) :: x,xx,norm,gam
        real(dp) :: sxpt,expnum,expden,tpg,sx,elwert
        real(dp) :: relfactor
        real(dp) :: gff
    
        rY = rho*ionic% Ye
        T8 = T*1.0e-8_dp
        T8_32 = T8**1.5

        ! test for overflow
        if (eta_e > log(huge(1.0_dp))) then
            xx  = eta_e
        else
            xx  = log(1.0+exp(eta_e))
        end if
  
        ! normalisation and degeneracy piece
        norm = 1.16_dp*8.02e3_dp*xx*T8_32/rY
  
        x = (1.0+xx)**twothird
        gam = sqrt(1.58e-3/T8)*ionic% Z;
  
      ! Elwert factor
        sxpt = sqrt(x+10.0)
        sx = sqrt(x)
        tpg = 2.0*pi*gam
        expnum = exp(-tpg/sxpt)
        expden = exp(-tpg/sx)
        elwert = (1.0 - expnum)/(1.0 - expden)

        ! relativistic piece
        relfactor = 1.0+(T8/7.7)**1.5
    
        ! gaunt factor and opacity
        gff = norm*elwert*relfactor
        kff = bfac*(rY/1.0e5_dp)/T8**3.5*gff*ionic% Z2/ionic% A
    end function freefree

    function Rosseland_kappa(rho,T,mu_e,ionic) result(kap)
        ! implements Rosseland mean for free-free and Thompson scattering 
        ! according to the fit of Potekhin and Yakovlev (2001, A&A 374: 213)
        !
        use nucchem_def, only: composition_info_type
        use constants_def
        real(dp), intent(in) :: rho,T,mu_e ! mu_e is in MeV
        type(composition_info_type), intent(in) :: ionic
        real(dp) :: kap
        real(dp) :: eta_e
        real(dp) :: T6, theta, kap_th, kap_ff, f, c7, T_Ry

        eta_e = mu_e*mev_to_ergs/boltzmann/T
        theta = T*boltzmann/Melectron/clight2
        kap_th = electron_scattering(eta_e,theta,ionic% Ye)
        kap_ff = freefree(rho,T,eta_e,ionic)
    
        T6 = T*1.0e-6_dp
        T_Ry = T6/0.15789_dp/ionic% Z2
    
        f = kap_ff/(kap_ff+kap_th)
        kap = (kap_th + kap_ff)*A(f,T_Ry)
    contains
        function A(f,T_Ry)
            real(dp), intent(in) :: f,T_Ry
            real(dp) :: A
            A = 1.0_dp + (1.097_dp+0.777_dp*T_Ry)/(1.0+0.536_dp*T_Ry)*f**0.617_dp*(1.0_dp-f)**0.77_dp
        end function A
    end function Rosseland_kappa

end module photon_conductivity
