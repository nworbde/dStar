module eval_conductivity
    use conductivity_def
    use superfluid_def, only: max_number_sf_types, neutron_1S0
    use superfluid_lib, only: sf_get_results    
    
    integer, parameter :: int_default_max_steps = 1000
    real(dp), parameter :: int_default_max_step_size = 0.0
    real(dp), parameter :: int_default_starting_step = 1.0e-3
    
    
contains
    
    subroutine conductivity(rho,T,chi,Gamma,eta,mu_e,ionic,kappa, &
        & which_ee,which_eQ,K_components)
        use constants_def
        use nucchem_def, only: composition_info_type
        real(dp), intent(in) :: rho, T, chi, Gamma, eta, mu_e   ! mu_e is in MeV
        type(composition_info_type), intent(in) :: ionic
        type(conductivity_components), intent(out) :: kappa
        integer, intent(in) :: which_ee, which_eQ
        integer :: ierr
        logical, dimension(num_conductivity_channels), intent(in) :: K_components
        real(dp) :: nn,nion, nu, nu_c, kappa_pre, ne, kF, xF, eF, Gamma_e, K_opacity
        real(dp) :: kappa_n_pre, nu_n, kfn, mnstar, tau, v, R, Tns
        real(dp) :: Tc(max_number_sf_types)
    
        ne = rho/amu*ionic%Ye
        kF = (threepisquare*ne)**onethird
        xF = hbar*kF/Melectron/clight
        eF = sqrt(1.0_dp+xF**2)
        Gamma_e = Gamma/ionic%Z53
        kappa_pre = onethird*pi**2*boltzmann**2*T*ne/Melectron/eF
    
        call clear_kappa
        nu_c = 0.0_dp
        nu = 0.0_dp
        K_opacity = 0.0_dp
        if (K_components(icond_ee)) then 
            if (which_ee == icond_sy06) then
                nu_c = ee_SY06(ne,T)
            else
                nu_c = ee_PCY(ne,T)
            end if
            nu = nu + nu_c
            kappa% ee = kappa_pre/nu_c
        end if
        if (K_components(icond_ei)) then
            nu_c = eion(kF,Gamma_e,eta,ionic%Ye,ionic%Z,ionic%Z2,ionic%Z53,ionic%A)
            nu = nu + nu_c
            if (nu_c > tiny(1.0_dp)) kappa% ei = kappa_pre/nu_c
        end if
        if (K_components(icond_eQ)) then
           if (which_eQ == icond_eQ_potekhin) then
               nu_c = eQ(kF,T,ionic%Ye,ionic%Z,ionic%Z2,ionic%A,ionic%Q)
            else
               nu_c = eQ_page(kF,T,ionic%Ye,ionic%Z,ionic%Z2,ionic%A,ionic%Q)
           end if
            nu = nu + nu_c
            if (ionic%Q > 1.0e-8_dp) then
                kappa% eQ = kappa_pre/nu_c
            else
                kappa% eQ = -1.0_dp
            end if
        end if
        kappa% electron_total = kappa_pre/nu
        if (K_components(icond_sf) .and. ionic% Yn > 0.0) then
            ! what if the neutrons aren't SF?
            nn = rho*ionic% Yn/(1.0_dp-chi)/Mneutron / density_n
            nion = (1.0_dp-ionic%Yn)*rho/Mneutron/ionic% A /density_n
            kappa% sf =  sPh(nn,nion,T,ionic)
        end if
        if (K_components(icond_kap)) then
            kappa% kap = Rosseland_kappa(rho,T,mu_e,ionic)
            K_opacity = 4.0_dp*onethird*arad*clight*T**3/rho/kappa% kap
        end if
        kappa% total = kappa% sf + K_opacity
        if (nu > 0.0_dp) kappa% total = kappa% total + kappa_pre/nu

        ! neutrons
        nu_n = 0.0_dp
        nu_c = 0.0_dp
        kappa_n_pre = 0.0_dp
    
        if (ionic% Yn > 0.0) then
            nn = rho/amu*ionic% Yn/(1.0-chi)
            kFn = (threepisquare*nn)**onethird
            nion = rho/amu*(1.0-ionic%Yn)/ionic%A   
            mnstar = Mneutron
            kappa_n_pre = onethird*pi**2*boltzmann**2*T*nn/mnstar    

            call sf_get_results(0.0_dp,kfn/cm_to_fm,Tc)
            ! Now set the superfluid reduction factor
            R = 1.0
            Tns = Tc(neutron_1S0)
            if (T < Tns) then
                tau = T/Tns
                v = sqrt(1.0-tau)*(1.456-0.157/sqrt(tau) + 1.764/tau)
                R = (0.4186+sqrt(1.007**2+(0.5010*v)**2))**2.5 * &
                &   exp(1.456-sqrt(1.456**2+v**2))
                if (R < 1.0E-10) R = 0.0
            end if

            ! note double multiplication by R
            if (K_components(icond_nQ)) then
                nu_c = n_imp(nn,nion,T,ionic,ierr)
                kappa% nQ = kappa_n_pre/nu_c*R
                nu_n = nu_n + nu_c
            end if
            if (K_components(icond_np)) then
                nu_c = n_phonon(nn,nion,T,ionic, ierr)
                kappa% np = kappa_n_pre/nu_c*R
                nu_n = nu_n + nu_c
            end if

            kappa % neutron_total = kappa_n_pre/nu_n*R

        end if     

        kappa % total = kappa% total + kappa % neutron_total

        contains
        subroutine clear_kappa()
            kappa% total = 0.0_dp
            kappa% ee  = 0.0_dp
            kappa% ei = 0.0_dp
            kappa% eQ = 0.0_dp
            kappa% sf = 0.0_dp
            kappa% kap = 0.0_dp
            kappa% nQ = 0.0_dp
            kappa% np = 0.0_dp
            kappa% electron_total = 0.0_dp
            kappa% neutron_total = 0.0_dp
        end subroutine clear_kappa
    end subroutine conductivity

    function ee_PCY(ne,T) result(nu)
        ! fmla. from Potekhin, Chabrier, and Yakovlev (1997)
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

        Il = (1.0/u)*(0.1587 - 0.02538/(1.0+0.0435*theta))*log(1.0+128.56/theta/(37.1 + theta*(10.83 + theta)))
        It = u3*(2.404/Ct + (Ct2-2.404/Ct)/(1.0+0.1*tu))*log(1.0 + Ct/tu/(At + tu))
        Ilt = u*(18.52*u2/Clt + (Clt2-18.2*u2/Clt)/(1.0+0.1558*theta**Blt)) &
                & *log(1.0 + Clt/(Alt*theta + 10.83*tu**2 + tu**(8.0/3.0)))

        I = Il + It + Ilt
        nu = nu_pre*ne*I/eF/T
    end function ee_SY06

    function eion(kF,Gamma_e,eta,Ye,Z,Z2,Z53,A) result(nu)
        ! implements fmla. of Baiko et al. (1999)
        use constants_def
        real(dp), intent(in) :: kF, Gamma_e, eta, Ye, Z, Z2, Z53, A
        real(dp) :: nu
        ! interface
        !   function eone(z)
        !       real(dp)(kind=8) :: z, eone
        !   end function eone
        ! end interface
        real(dp), parameter :: um1 = 2.8, um2 = 12.973, onesixth = 1.0/6.0, fourthird = 4.0*onethird
        real(dp) :: electroncompton
        real(dp) :: aB
        real(dp) :: eifac
        real(dp) :: kF2,x,eF,Gamma,v,v2,kTF2,eta02,aei,qD2,beta,qi2,qs2
        real(dp) :: Gs,Gk0,Gk1,GkZ,Gk,a0,D1,D,s,w1,w,sw,fac
        real(dp) :: L1,L2,Lei
    
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
        fac = exp(sw)*(eone(sw)-eone(sw+w))
    
        L1 = 0.5*(log(1.0+1.0/s) + (1.0-exp(-w))*s/(s+1.0) - fac*(1.0+sw))
        L2 = 0.5*(1.0-(1.0-exp(-w))/w + s*(s/(1.0+s)*(1.0-exp(-w))  &
                    & - 2.0*log(1.0+1.0/s)  +fac*(2.0+sw)))
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

    function sPh(nn, nion, temperature, ionic)
        ! neutron superfluid phonon conductivity, following Aguilera et al. (2009), PRL
        !
        use nucchem_def, only: composition_info_type
        use constants_def
        real(dp), intent(in) :: nn, nion, temperature
        type(composition_info_type), intent(in) :: ionic
        real(dp) :: sPh
        real(dp) :: K_n
        real(dp) :: etwo
        real(dp) :: Mu_n
        real(dp), parameter :: anl = 10.0
        real(dp) :: T,ai,kFn,kFe,vs,Cv,omega_p,TUm,anlt,gmix,cs,qTFe,alpha,tau_lph,wt
        real(dp) :: Llphn,Llphu,Llph,Lsph,fu,omega

        K_n = boltzmann*clight*fm_to_cm/(hbarc_n*fm_to_cm)**3
        etwo = finestructure*hbarc_n
        Mu_n = amu*clight2*ergs_to_mev
    
        ! convert to nuclear units
        T = temperature*boltzmann*ergs_to_mev
        ai = (3.0/fourpi/nion)**onethird
        kFn = (threepisquare*nn)**onethird
        kFe = (threepisquare*nion*ionic%Z)**onethird
        ! phonon speed
        vs = hbarc_n*kFn/Mn_n/sqrt(3.0)
        ! specific heat, phonon gas
        Cv = 2.0*pi**2*T**3/15.0/vs**3
        ! ion plasma freq
        omega_p = sqrt(fourpi*etwo*ionic%Z2*nion/ionic%A/Mu_n)
        ! Thomas-Fermi wavevector
        qTFe = sqrt(4.0*finestructure/pi)*kFe
        ! sound speed
        cs = omega_p/qTFe
        alpha = cs/vs
        ! mixing parameter
        anlt = anl/(1.0+0.4*kFn*anl + 0.5*2.0*anl*kFn**2)
        gmix = 2.0*anlt*hbarc_n*sqrt(nion*kFn/ionic%A/Mu_n**2)
        ! Umklapp temperature
        TUm = onethird*etwo*omega_p*(ionic% Z)**onethird
        ! thermal phonon frequency
        omega = 3.0*T/hbarc_n
        ! interpolation fmla.
        fu = exp(-TUm/T)
        ! mean-free path, normal; fu interpolates to Umklapp scattering
        Llphn = 2.0/pi/omega
        Llphu = 100.0*ai
        Llph = 1.0/(fu/Llphu + (1.0-fu)/Llphn)
        tau_lph = Llph/cs
        wt = tau_lph*omega
        Lsph = Llph*(vs/gmix)**2*(1.0+(1.0-alpha**2)**2 * wt**2)/alpha/wt**2
        sPh = onethird*Cv*vs*Lsph * K_n
    end function sPh

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
        kTh = 8.0*onethird*pi*avogadro*(electroncharge**2/Melectron/clight2)**2*Ye/Gbar
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
        ! implements Rosseland mean for free-free and Thompson scattering according to the fit of
        ! Potekhin and Chabrier (2001, A&A)
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


    function n_imp(nn, nion, temperature, ionic, ierr) result(nu)
        ! neutron-impurity scattering, following S. Reddy's notes
        use, intrinsic :: iso_fortran_env, only: error_unit
        use nucchem_def, only: composition_info_type
        use constants_def
        use num_lib

        real(dp), intent(in) :: nn, nion, temperature
        type(composition_info_type), intent(in) :: ionic
        real(dp) :: ne      ! number density of electrons
        real(dp) :: kFn     ! neutron Fermi wavevector 
        real(dp) :: EFn     ! neutron Fermi energy
        real(dp) :: nu      ! neutron-impurity scattering frequency
        real(dp) :: V       ! effective neutron-impurity potential 
        real(dp) :: Lambda_n_imp    ! Coulomb logarithm (Potekhin et al. 1999)
        real(dp) :: fac, sf_frac, Tc, xFn
        real(dp) :: gamma, vFn, v2, s, L1
        real(dp) :: pfn
        real(dp) :: R_a, mnstar
        real(dp) :: isospin, n_0, n_2, n_l, n_in    
        real(dp), dimension(:), pointer :: y
        integer, intent(out) :: ierr
        integer ,dimension(:), pointer :: iwork => null()
        real(dp), dimension(:), pointer :: work => null()
        real(dp), dimension(1) :: rtol, atol
        integer, dimension(:), pointer :: ipar => null()
        real(dp), dimension(:), pointer :: rpar => null()
        real(dp) :: h, kn_end, kn_start, kn, kfi
        integer :: liwork, lwork, itol, lipar, lrpar
        integer :: n, idid, lout, iout, npts, num_integration_variables, num_ipar, num_rpar

        allocate(y(1))
        num_ipar = 0
        num_rpar = 6    
        kFn = (threepisquare*nn)**onethird
        kfi = (threepisquare*nn/ionic%Yn*(1.0-ionic%Yn)/ionic%A)**onethird
        num_integration_variables = 1
        call dop853_work_sizes( &
        &   num_integration_variables,num_integration_variables,liwork,lwork)
        allocate(iwork(liwork), work(lwork),ipar(num_ipar), rpar(num_rpar))
        iwork = 0
        iwork(5) = num_integration_variables
        work = 0.0
        
        itol = 0
        rtol = 1.0e-4
        atol = 1.0e-5
        iout = 0    ! want dense output
        lout = error_unit
        lipar = num_ipar
        lrpar = num_rpar
        n = num_integration_variables
                    
        !calculate n_in and R_a using Sly4 model for now
        isospin = 1.0 - 2.0*ionic%Z/ionic%A
        n_0 = 0.1740/(fm_to_cm)**3
        n_2 = -0.0157/(fm_to_cm)**3
        n_l = n_0+n_2*isospin**2
        n_in = n_l/2.0*(1.0+0.9208*isospin)
        R_a = (3.0*(ionic%A-ionic%Z)/4.0/pi/n_in)**onethird    
    
        y(1) = 0.0
        kn_start = 1.0E-10 
        kn_end = 2.0*kFn*R_a
        h = -0.1

        rpar(1) = temperature
        rpar(2) = ionic% Z
        rpar(3) = ionic% A
        rpar(4) = ionic% Yn
        rpar(5) = kfn
        rpar(6) = R_a

        call dop853(n,structure_factor_impurity,kn_start,y,kn_end,h, &
        &   int_default_max_step_size,int_default_max_steps, &
        & rtol,atol,itol, solout, iout, work, lwork, iwork, liwork,  &
        &   num_rpar, rpar, num_ipar, ipar, lout, idid)

        deallocate(work, iwork)
    
        Lambda_n_imp = y(1)
        fac = 4.0/27.0/pi/hbar**3/clight**2
        V = neutron_potential(nn,n_in,ionic%Z,ionic%A)    
        kFn = (threepisquare*nn)**onethird
        xFn = hbar*kFn/Mneutron/clight
        mnstar = Mneutron   
        Efn = 15.1*mev_to_ergs*(ionic% Yn)**(2.0/3.0)* &
        &   (nn*amu/ionic% Yn/1.0E14)**(2.0/3.0)* &
        &   (2.0*Mneutron/mnstar) + Mneutron*clight**2
        
        nu = fac*EFn*(nion/nn)*(V*R_a)**2*ionic%Q/(ionic%Z)**2*Lambda_n_imp
    end function n_imp

    function n_phonon(nn, nion, temperature, ionic, ierr) result(nu)
        ! neutron-phonon scattering, following S. Reddy's notes
        !
        use, intrinsic :: iso_fortran_env, only: error_unit
        use nucchem_def, only: composition_info_type
        use constants_def
        use num_lib    
    
        real(dp), intent(in) :: nn, nion, temperature
        type(composition_info_type), intent(in) :: ionic
        real(dp) :: ne      ! electron number density
        real(dp) :: kFn     ! neutron Fermi wavevector 
        real(dp) :: EFn     ! neutron Fermi energy
        real(dp) :: nu      ! neutron-impurity scattering frequency
        real(dp) :: V       ! effective neutron-impurity potential 
        real(dp) :: Lambda_n_phonon     ! Coulomb logarithm (Potekhin et al. 1999)
        real(dp) :: fac, sf_frac, Tc, xFn
        real(dp) :: kfe, xr_e, eta, R_a, mnstar
        real(dp) :: qD, qi, ktf, beta, qs, wf, sf, L2, eta_n_0, G_k, Tp, eta_n
        real(dp) :: isospin, n_0, n_2, n_l, n_in
        real(dp), dimension(:), pointer :: y
        integer, intent(out) :: ierr
        integer ,dimension(:), pointer :: iwork => null()
        real(dp), dimension(:), pointer :: work => null()
        real(dp), dimension(1) :: rtol, atol
        integer, dimension(:), pointer :: ipar => null()
        real(dp), dimension(:), pointer :: rpar => null()
        real(dp) :: h, kn_end, kn_start, kn, kfi
        integer :: liwork, lwork, itol, lipar, lrpar
        integer :: n, idid, lout, iout, npts, num_integration_variables, num_ipar, num_rpar

        allocate(y(1))
        num_ipar = 0
        num_rpar = 6    
        kFn = (threepisquare*nn)**onethird
        kfi = (threepisquare*nn/ionic%Yn*(1.0-ionic%Yn)/ionic%A)**onethird
        num_integration_variables = 1
        call dop853_work_sizes( &
        &   num_integration_variables,num_integration_variables,liwork,lwork)
        allocate(iwork(liwork), work(lwork),ipar(num_ipar), rpar(num_rpar))
        iwork = 0
        iwork(5) = num_integration_variables
        work = 0.0
        
        itol = 0
        rtol = 1.0e-4
        atol = 1.0e-5
        iout = 0    ! want dense output
        lout = error_unit
        lipar = num_ipar
        lrpar = num_rpar
        n = num_integration_variables
                  
        y(1) = 0.0

        !calculate n_in and R_a using Sly4 model for now
        isospin = 1.0 - 2.0*ionic%Z/ionic%A
        n_0 = 0.1740/(fm_to_cm)**3
        n_2 = -0.0157/(fm_to_cm)**3
        n_l = n_0+n_2*isospin**2
        n_in = n_l/2.0*(1.0+0.9208*isospin)
        R_a = (3.0*(ionic%A-ionic%Z)/4.0/pi/n_in)**onethird
    
        kn_start = 1.0E-10 
        kn_end = 2.0*kFn*R_a
        h = -0.1
        rpar(1) = temperature
        rpar(2) = ionic% Z
        rpar(3) = ionic% A
        rpar(4) = ionic% Yn
        rpar(5) = kfn
        rpar(6) = R_a

        call dop853(n,structure_factor_phonon,kn_start,y,kn_end,h, &
        &   int_default_max_step_size,int_default_max_steps, &
        &   rtol,atol,itol, solout, iout, work, lwork, iwork, liwork,  &
        &   num_rpar, rpar, num_ipar, ipar, lout, idid)

        deallocate(work, iwork)
    
        Lambda_n_phonon = y(1)
        fac = 4.0/27.0/pi/hbar**3/clight**2
        V = neutron_potential(nn,n_in,ionic%Z,ionic%A) ! V in [ergs]
        kFn = (threepisquare*nn)**onethird
        xFn = hbar*kFn/Mneutron/clight
        mnstar = Mneutron
        Efn = 15.1*mev_to_ergs* &
        &   (ionic% Yn)**(2.0/3.0)*(nn*amu/ionic%Yn/1.0E14)**(2.0/3.0)* &
        &   (2.0*Mneutron/mnstar) + Mneutron*clight**2 

        nu = fac*EFn*(nion/nn)*(V*R_a)**2*Lambda_n_phonon
    end function n_phonon

    function neutron_potential(nn,n_in,Z,A) result(V)
        use constants_def
        real(dp), intent(in) :: nn, n_in, Z, A
        real(dp) :: V
        V = hbar**2*(threepisquare*n_in)**(2.0/3.0)/2.0/Mneutron &
        &   *(1.0 - (nn/n_in)**(2.0/3.0))
    end function neutron_potential


    subroutine structure_factor_phonon(n,x,h,y,dy,lrpar,rpar,lipar,ipar,ierr)
        use constants_def

        integer, intent(in) :: n, lrpar, lipar
        real(dp), intent(in) ::  x,h
        real(dp), intent(inout) :: y(:)
        real(dp), intent(out) :: dy(:)
        real(dp), intent(inout), pointer :: rpar(:)
        integer, intent(inout), pointer :: ipar(:)
        integer, intent(out) :: ierr
        real(dp) :: Yn, nn, nion, kfn
        real(dp) :: mion, omegaP, eta, Gamma, T, aion, A, Z, q, formfac
        real(dp) :: a1, um1, um2, w, w1, ssig, dskappa2, dskappa4, dskappa
        real(dp) :: R_a

        ierr = 0
        T = rpar(1)
        Z = rpar(2) 
        A = rpar(3)
        Yn = rpar(4)
        kfn = rpar(5)
        R_a = rpar(6)
        
        nn = kfn**3/threepisquare
        nion = nn/Yn*(1.0-Yn)/A 
        aion = (3.0/(4.0*pi*nion))**onethird
        mion = (A-Z)*Mneutron+Z*Mproton
        omegaP=(hbar*clight*4.*pi*Z**2*finestructure*nion/mion)**(1./2.)
        eta = boltzmann*T/(hbar*omegaP)
        Gamma = hbar*clight*Z**2*finestructure/(boltzmann*T*aion)
        a1= (x/R_a*aion)**2/(3.*Gamma*eta)
        um1=2.8
        um2=13.0 
        w = a1*(um1*exp(-9.1*eta)+2.*eta*um2)/4.
        w1= a1*um2*eta**2/(2.*(eta**2+(um2/117.)**2)**(1./2.))
        ssig=exp(-2.*(w-w1))*(1.0-exp(-2.*w1))
        dskappa2 = a1*91.*eta**2*exp(-2.*w)/(1.+111.4*eta**2)**2 
        dskappa4 = a1*0.101*eta**4/((0.06408+eta**2)*(0.001377+eta**2)**1.5)
        dskappa = dskappa2+dskappa4        
        formfac = abs(3.0*(dsin(x)-(x)*dcos(x))/(x)**3)
        
        dy(1) = x**3*formfac**2*(ssig+(3.*(kfn*R_a/x)**2-0.5)*dskappa)
        
    end subroutine structure_factor_phonon

    subroutine structure_factor_impurity(n,x,h,y,dy,lrpar,rpar,lipar,ipar,ierr)
        use constants_def

        integer, intent(in) :: n, lrpar, lipar
        real(dp), intent(in) ::  x,h
        real(dp), intent(inout) :: y(:)
        real(dp), intent(out) :: dy(:)
        real(dp), intent(inout), pointer :: rpar(:)
        integer, intent(inout), pointer :: ipar(:)
        integer, intent(out) :: ierr
        real(dp) :: Yn, nn, nion, kfn
        real(dp) :: mion, omegaP, eta, Gamma, T, aion, A, Z, q, formfac
        real(dp) :: a1, um1, um2, w, w1, ssig, dskappa2, dskappa4, dskappa
        real(dp) :: R_a

        ierr = 0

        T = rpar(1)
        Z = rpar(2) 
        A = rpar(3)
        Yn = rpar(4)
        kfn = rpar(5)
        R_a = rpar(6)

        nn = kfn**3/threepisquare
        nion = nn/Yn*(1.0-Yn)/A         
        aion = (3.0/(4.0*pi*nion))**onethird
        mion = (A-Z)*Mneutron+Z*Mproton
        omegaP=(hbar*clight*4.*pi*Z**2*finestructure*nion/mion)**(1./2.)
        eta = boltzmann*T/(hbar*omegaP)
        Gamma = hbar*clight*Z**2*finestructure/(boltzmann*T*aion)
        a1= (x/R_a*aion)**2/(3.*Gamma*eta)
        um1=2.8
        um2=13.0 
        w = a1*(um1*exp(-9.1*eta)+2.*eta*um2)/4.
        w1= a1*um2*eta**2/(2.*(eta**2+(um2/117.)**2)**(1./2.))
        ssig=exp(-2.*(w-w1))*(1.0-exp(-2.*w1))
        dskappa2 = a1*91.*eta**2*exp(-2.*w)/(1.+111.4*eta**2)**2 
        dskappa4 = a1*0.101*eta**4/((0.06408+eta**2)*(0.001377+eta**2)**1.5)
        dskappa = dskappa2+dskappa4        
        formfac = abs(3.0*(dsin(x)-(x)*dcos(x))/(x)**3)
        
        dy(1) = x**3*formfac**2
        
    end subroutine structure_factor_impurity

    subroutine solout(nr, xold, x, n, y, rwork_y, iwork_y, interp_y, lrpar, rpar, lipar, ipar, irtrn)
        use dStar_crust_lib
        
        integer, intent(in) :: nr, n, lrpar, lipar
        real(dp), intent(in) :: xold, x
        real(dp), intent(inout) :: y(:)
        real(dp), intent(inout), target :: rwork_y(*)
        integer, intent(inout), target :: iwork_y(*)
        integer, intent(inout), pointer :: ipar(:) ! (lipar)
        real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
        integer, intent(out) :: irtrn
        interface
            double precision function interp_y(i, s, rwork_y, iwork_y, ierr)
            integer, intent(in) :: i ! result is interpolated approximation of y(i) at x=s.
            double precision, intent(in) :: s ! interpolation x value (between xold and x).
            double precision, intent(inout), target :: rwork_y(*)
            integer, intent(inout), target :: iwork_y(*)
            integer, intent(out) :: ierr
         end function interp_y
        end interface ! 
        integer :: ierr, i

   !     if (ierr /= 0) irtrn = -1
    end subroutine solout


end module eval_conductivity
