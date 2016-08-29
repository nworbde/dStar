module eval_conductivity
    use conductivity_def
    
contains

subroutine conductivity(rho,T,chi,Gamma,eta,mu_e,ionic,kappa,which_ee,which_eQ,K_components)
    use constants_def
    use nucchem_def, only: composition_info_type
    real(dp), intent(in) :: rho, T, chi, Gamma, eta, mu_e   ! mu_e is in MeV
    type(composition_info_type), intent(in) :: ionic
    type(conductivity_components), intent(out) :: kappa
    integer, intent(in) :: which_ee, which_eQ
    logical, dimension(num_conductivity_channels), intent(in) :: K_components
<<<<<<< HEAD
    real(dp) :: nn,nion, nu, nu_c, kappa_pre, ne, kF, xF, eF, Gamma_e
    real(dp) :: nu_n, kappa_n_pre, kFn, xFn, eFn
=======
    real(dp) :: nn,nion, nu, nu_c, kappa_pre, ne, kF, xF, eF, Gamma_e, K_opacity
>>>>>>> nworbde/master
    
    ne = rho/amu*ionic%Ye
    kF = (threepisquare*ne)**onethird
    xF = hbar*kF/Melectron/clight
    eF = sqrt(1.0+xF**2)
    Gamma_e = Gamma/ionic%Z53
    kappa_pre = onethird*pi**2*boltzmann**2*T*ne/Melectron/eF
    kappa_n_pre = 0.0
        
    ! electrons
    call clear_kappa
    nu_c = 0.0
    nu = 0.0
    K_opacity = 0.0
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
        if (ionic%Q > 1.0e-8) then
            kappa% eQ = kappa_pre/nu_c
        else
            kappa% eQ = -1.0
        end if
    end if
    kappa% total = kappa_pre/nu 
    
    ! neutrons
    nu_n = 0.0  
    if (ionic% Yn > 0.0) then
    nn = rho*ionic% Yn/(1.0-chi)/Mneutron / density_n
    kFn = (threepisquare*nn)**onethird
    xFn = hbar*kFn/Mneutron/clight
    eFn = sqrt(1.0+xFn**2)
    nion = (1.0-ionic%Yn)*rho/Mneutron/ionic% A /density_n    
    kappa_n_pre = onethird*pi**2*boltzmann**2*T*nn/Mneutron/eFn      
    if (K_components(icond_sf)) then
        kappa% sf =  sPh(nn,nion,T,ionic)
    end if
<<<<<<< HEAD
	if (K_components(icond_nQ)) then
    	nu_c = n_imp(nn,nion,T,ionic)
    	kappa% nQ = kappa_n_pre/nu_c
        nu_n = nu_n + nu_c
    end if    
 	if (K_components(icond_np)) then
        nu_c = n_phonon(nn,nion,T,ionic)
        kappa% np = kappa_n_pre/nu_c
        nu_n = nu_n + nu_c
    end if      
		kappa% total = kappa% total + kappa_n_pre/nu_n + kappa% sf
	end if
=======
    if (K_components(icond_kap)) then
        kappa% kap = Rosseland_kappa(rho,T,mu_e,ionic)
        K_opacity = 4.0*onethird*arad*clight*T**3/rho/kappa% kap
    end if
    kappa% total = kappa% sf + K_opacity
    if (nu > 0.0) kappa% total = kappa% total + kappa_pre/nu
>>>>>>> nworbde/master
    
    contains
    subroutine clear_kappa()
        kappa% total = 0.0
        kappa% ee  = 0.0
        kappa% ei = 0.0
        kappa% eQ = 0.0
<<<<<<< HEAD
        kappa% sf = 0.0     
        kappa% nQ = 0.0
        kappa% np = 0.0
=======
        kappa% sf = 0.0    
        kappa% kap = 0.0 
>>>>>>> nworbde/master
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

<<<<<<< HEAD
function n_imp(nn, nion, temperature, ionic) result(nu)
    ! neutron-impurity scattering, following S. Reddy's notes
    !
    use nucchem_def, only: composition_info_type
    use constants_def
    real(dp), intent(in) :: nn, nion, temperature
    type(composition_info_type), intent(in) :: ionic
    real(dp) :: ne		! number density of electrons
    real(dp) :: kFn 	! neutron Fermi wavevector 
    real(dp) :: EFn 	! neutron Fermi energy
    real(dp) :: nu		! neutron-impurity scattering frequency
    real(dp) :: V_0 	! effective neutron-impurity potential at V(R)
    real(dp) :: R_a		! radius of scattering center
    real(dp) :: Lambda_n_imp 	! Coulomb logarithm (Potekhin et al. 1999)
	real(dp) :: fac, sf_frac, Tc

	V_0 = 20.0*mev_to_ergs
	R_a = 10.0*fm_to_cm
	Tc = 1.0E8	
	Lambda_n_imp = 1.0
    fac = 4.0/27.0/pi/hbar**2/clight
    sf_frac = exp(sqrt(1.0-temperature/Tc))
   
    kFn = (3.0*pi*nn)**onethird
    EFn = hbar*clight*kFn
    ne = ionic%Z*nion 
    nu = fac*EFn*(ne/nn)*(V_0*R_a)**2*Lambda_n_imp*ionic%Q/(ionic%Z)**3
end function n_imp

function n_phonon(nn, nion, temperature, ionic) result(nu)
    ! neutron-phonon scattering, following S. Reddy's notes
    !
    use nucchem_def, only: composition_info_type
    use constants_def
    real(dp), intent(in) :: nn, nion, temperature
    type(composition_info_type), intent(in) :: ionic
    real(dp) :: ne		! electron number density
    real(dp) :: kFn 	! neutron Fermi wavevector 
    real(dp) :: EFn 	! neutron Fermi energy
    real(dp) :: nu		! neutron-impurity scattering frequency
    real(dp) :: V_0 	! effective neutron-impurity potential at V(R)
    real(dp) :: R_a		! radius of scattering center
    real(dp) :: Lambda_n_phonon 	! Coulomb logarithm (Potekhin et al. 1999)
	real(dp) :: fac, sf_frac, Tc

	V_0 = 20.0*mev_to_ergs
	R_a = 10.0*fm_to_cm
	Tc = 1.0E8
	Lambda_n_phonon = 1.0
    fac = 4.0/27.0/pi/hbar**2/clight
    ! suppression near critical temperature
    sf_frac = exp(sqrt(1.0-temperature/Tc))
   
    kFn = (3.0*pi*nn)**onethird
    EFn = hbar*clight*kFn
    ne = ionic%Z*nion 
    nu = fac*EFn*(nion/nn)*(V_0*R_a)**2*Lambda_n_phonon
end function n_phonon
=======
function electron_scattering(eta_e,theta,Ye) result(kTh)
    use constants_def
    real(dp), intent(in) :: eta_e   ! electron chemical pot./kT
    real(dp), intent(in) :: theta ! T/mc**2
    real(dp), intent(in) :: Ye
    real(dp) :: kTh
    real(dp) :: xi, xi2, theta2
    real(dp) :: t1, t2, t3, Gbar
    
    xi	= exp(0.8168*eta_e - 0.05522*eta_e**2)
    xi2	= xi**2
    theta2 = theta**2
    t1	= 1.129 + 0.2965*xi - 0.005594*xi2
    t2	= 11.47 + 0.3570*xi + 0.1078*xi2
    t3	= -3.249 + 0.1678*xi - 0.04706*xi2
  
    Gbar	=  t1 + t2*theta + t3*theta2
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
>>>>>>> nworbde/master

end module eval_conductivity
