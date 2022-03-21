module electron_conductivity
    use conductivity_def
    implicit none
    
contains
    
    function ee_PCY(ne,T) result(nu)
        ! fmla. from Potekhin, Chabrier, and Yakovlev (1997)
        use math_lib
        use constants_def
        real(dp), intent(in) :: ne, T
        real(dp) :: nu
        real(dp) :: yfac
        real(dp) :: mec2
        real(dp) :: eefac
        real(dp) :: befac
        real(dp) :: onesixth = 1.0/6.0
        real(dp) :: plasma_theta, kF
        real(dp) :: x, eF, beta, b2, be3, y, lJ, llJ, J0,J1,J2,J

        yfac = sqrt(4.0*finestructure/pi)
        mec2 = Melectron*clight**2
        eefac = 1.5*finestructure**2*mec2/pi**3/hbar
        befac = (finestructure/pi)**1.5
    
        kF = (3.0*pi**2*ne)**onethird
        x = kF*hbar/Melectron/clight
        eF = sqrt(1.0+x**2)
        beta = x/eF
        b2 = beta**2
        be3 = befac/beta**1.5
        plasma_theta = boltzmann*T/mec2
        y = yfac*x**1.5/plasma_theta/sqrt(ef)
        lJ = 1.0 + (2.810-0.810*b2)/y
        llj = log(lJ)
        J0 = onethird/(1.0+0.07414*y)**3
        J1 = J0*llJ + onesixth*pi**5*y/(13.91+y)**4
        J2 = 1.0+0.4*(3.0+1.0/x**2)/x**2
        J = y**3 * J1*J2
    
        nu = eefac*plasma_theta**2*J/eF/be3
    end function ee_PCY

    function ee_SY06(ne,T) result(nu)
        ! uses fmla of Shternin & Yakovlev (2006). Checked against routine courtesy of D. Page.
        !
        use math_lib
        use constants_def
        real(dp), intent(in) :: ne, T
        real(dp) :: nu
        real(dp) :: nu_pre
        real(dp) :: Tp_pre
        real(dp) :: Il, It, Ilt, I
        real(dp) :: At, Ct, Ct1, Ct2, Alt, Blt, Clt, Clt1, Clt2
        real(dp) :: u2, u3, u4, tu
        real(dp) :: kF, Tp, theta, xF, eF, u, mstar

        nu_pre = 36.0*finestructure**2*hbar**2*clight/pi/boltzmann/Melectron
        Tp_pre = hbar*sqrt(4.0*pi*finestructure*hbar*clight/Melectron)/boltzmann

        kF = (threepisquare*ne)**onethird
        xF = hbar*kF/Melectron/clight
        eF = sqrt(1.0+xF**2)
        u = xF/eF
        Tp = Tp_pre*sqrt(ne/eF)
        theta = sqrt(3.0)*Tp/T

        u2 = u**2
        u3 = u**3
        u4 = u**4
        tu = theta*u
        At = 20.0+450.0*u3
        Ct1 = 0.05067+0.03216*u2
        Ct2 = 0.0254+0.04127*u4
        Ct = At*exp(Ct1/Ct2)

        Alt = 12.2 + 25.2*u3
        Blt = 1.0-0.75*u
        Clt1 = 0.123636+0.016234*u2
        Clt2 = 0.0762+0.05714*u4
        Clt = Alt*exp(Clt1/Clt2)

        Il = (1.0/u)*(0.1587 - 0.02538/(1.0+0.0435*theta)) &
        &       *log(1.0+128.56/theta/(37.1 + theta*(10.83 + theta)))
        It = u3*(2.404/Ct + (Ct2-2.404/Ct)/(1.0+0.1*tu)) &
        &       *log(1.0 + Ct/tu/(At + tu))
        Ilt = u*(18.52*u2/Clt + (Clt2-18.2*u2/Clt)/(1.0+0.1558*theta**Blt)) &
        &       *log(1.0 + Clt/(Alt*theta + 10.83*tu**2 + tu**(8.0/3.0)))

        I = Il + It + Ilt
        nu = nu_pre*ne*I/eF/T
    end function ee_SY06

    function eion(kF,Gamma_e,eta,Ye,Z,Z2,Z53,A) result(nu)
        ! implements fmla. of Baiko et al. (1999)
        use math_lib
        use constants_def
        real(dp), intent(in) :: kF, Gamma_e, eta, Ye, Z, Z2, Z53, A
        real(dp) :: nu
        ! interface
        !   function eone(z)
        !       real(dp)(kind=8) :: z, eone
        !   end function eone
        ! end interface
        real(dp), parameter :: um1 = 2.8, um2 = 12.973
        real(dp), parameter :: onesixth = 1.0/6.0, fourthird = 4.0*onethird
        real(dp) :: sw_max
        real(dp) :: electroncompton
        real(dp) :: aB
        real(dp) :: eifac
        real(dp) :: kF2,x,eF,Gamma,v,v2,kTF2,eta02,aei,qD2,beta,qi2,qs2
        real(dp) :: Gs,Gk0,Gk1,GkZ,Gk,a0,D1,D,s,w1,w,sw,fac
        real(dp) :: L1,L2,Lei
    
        sw_max = log(huge(1.0_dp))
        electroncompton = hbar/Melectron/clight
        aB = electroncompton/finestructure
        eifac = fourthird*clight*finestructure**2/electroncompton/pi
    
        kF2 = kF**2
        x = electroncompton*kF
        eF = sqrt(1.0+x**2)
        v = x/eF
        v2 = v**2
        Gamma = Gamma_e*Z53
        s = finestructure/pi/v
        kTF2 = 4.0*s*kF2
        eta02 = (0.19/Z**onesixth)**2
        aei = (4.0/9.0/pi)**onethird * kF
        qD2 = 3.0*Gamma_e*aei**2 * Z2/Z
        beta = pi*finestructure*Z*v

        qi2 = qD2*(1.0+0.06*Gamma)*exp(-sqrt(Gamma))
        qs2 = (qi2+kTF2)*exp(-beta)
    
        Gs = eta/sqrt(eta**2+eta02)*(1.0+0.122*beta**2)
        Gk0 = 1.0/(eta**2+0.0081)**1.5
        Gk1 = 1.0 + beta*v**3
        GkZ = 1.0-1.0/Z
        Gk = Gs + 0.0105*GkZ*Gk1*Gk0*eta
        
        a0 = 4.0*kF2/qD2/eta
        D1 = exp(-9.1*eta)
        D = exp(-a0*um1*D1*0.25)
    
        w1 = 1.0+onethird*beta
        w = 4.0*um2*kF2/qD2*w1
    
        sw = s*w
        if (sw < sw_max) then
            fac = exp(sw)*(eone(sw)-eone(sw+w))
            L1 = 0.5*(log(1.0+1.0/s) + (1.0-exp(-w))*s/(s+1.0) - fac*(1.0+sw))
            L2 = 0.5*(1.0-(1.0-exp(-w))/w + s*(s/(1.0+s)*(1.0-exp(-w))  &
            & - 2.0*log(1.0+1.0/s)  +fac*(2.0+sw)))
        else
            L1 = 0.5*(log((s+1.0)/s) - 1.0/(s+1.0))
            L2 = (2.0*s+1.0)/(2.0*s+2.0) - s*log((s+1.0)/s)
        end if
        Lei = (Z2/A)*(L1-v2*L2)*Gk*D
        nu = eifac*eF*Lei*A/Z
    
        contains
        function eone(z)
            real(dp), intent(in) :: z
            real(dp) :: eone
        
            real(dp), parameter :: a0 = 0.1397, a1 = 0.5772, a2 = 2.2757
            real(dp) :: z2, z3, z4
            z2 = z*z; z3 = z2*z; z4 = z2*z2
            eone = exp(-z4/(z3+a0))*(log(1.0+1.0/z)-a1/(1.0+a2*z2))
        end function eone
    end function eion

    function eQ(kF,T,Ye,Z,Z2,A,Q) result (nu)
        ! impurity scattering: Iton (1994), also Potekhin, priv. communication
        use math_lib
        use constants_def
        real(dp), intent(in) :: kF,T,Ye,Z,Z2,A,Q
        real(dp) :: nu
        real(dp) :: dZ2,kF2,x,gamma,v2,v,qs2,s,L1,L
        real(dp) :: electron_compton, mec2
        real(dp) :: fac

        electron_compton = hbar/Melectron/clight
        mec2 = Melectron*clight**2
        fac = 4.0*onethird*clight*finestructure**2/electron_compton/pi
    
        dZ2 = Q/A
        kF2 = kF**2
        x = electron_compton*kF
        gamma = sqrt(1.0+x**2)
        v = x/gamma
        v2 = v**2
        s = finestructure/pi/v
        qs2 = 4.0*kF2*s
    
        L1 = 1.0+4.0*v2*s
        L = 0.5*(L1*log(1.0+1.0/s) - v2 -1.0/(s+1.0) - v2*s/(1.0+s))
        nu = fac*gamma*dZ2*L*A/Z
    end function eQ

    function eQ_page(kF,T,Ye,Z,Z2,A,Q) result(nu)
        use math_lib
        use constants_def
        real(dp), intent(in) :: kF,T,Ye,Z,Z2,A,Q
        real(dp) :: nu
        real(dp) :: dZ2,kF2,x,v,qs,L,gamma
        real(dp) :: electron_compton, mec2, fac
        
        electron_compton = hbar/Melectron/clight
        mec2 = Melectron*clight**2
        fac = 4.0*onethird*clight*finestructure**2/electron_compton/pi
        
        dZ2 = Q/A
        kF2 = kF**2
        x = electron_compton*kF
        gamma = sqrt(1.0+x**2)
        v = x/gamma
        qs = 0.048196/sqrt(v)
        L = log(1.0/qs)-0.5*(1+v**2)
        nu = fac*gamma*dZ2*L*A/Z
    end function eQ_page

end module electron_conductivity
