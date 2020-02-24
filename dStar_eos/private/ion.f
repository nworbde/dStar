! implementation of Chabrier & Potekhin (1998, PRE 58, 4941) fmla. for liquid;
! with Chabrier & Potekhin (2000) corrections
! Baiko, Potekhin, & Yakovlev (2001, PRE 64, 057402) fmla. for lattice (bcc)
! with Potekhin & Chabrier (2010) corrections
!
! Much of the code here are customized versions of the publically available routines 
! made available by A.Y. Potekhin at http://www.ioffe.ru/astro/EIP/
!

module ion_eos

    ! for returning ion lattice results
    integer, parameter :: number_lattice_products = 3
    integer, parameter :: total = 1
    integer, parameter :: harmonic = 2
    integer, parameter :: anharmonic = 3

contains

subroutine ion_mixture(rq,rs,Gamma_e,ionic,ncharged,charged_ids,Yion,Gamma,ionQ,phase,f,u,p,s,cv,dpr,dpT,err)
    use, intrinsic :: iso_fortran_env, only: error_unit
    use constants_def, only: dp, amu, Melectron
    use dStar_eos_def
    use dStar_eos_private_def
    use nucchem_def
    
    type(dStar_eos_general_info), pointer :: rq
    real(dp), intent(in) :: rs,Gamma_e
    type(composition_info_type), intent(in) :: ionic
    integer, intent(in) :: ncharged
    integer, dimension(ncharged), intent(in) :: charged_ids
        ! ids of the charged species
    real(dp), dimension(ncharged), intent(in) :: Yion
        ! renormalized abunances of charged species Yion = Y/(1-Yn)
    real(dp), intent(out) :: Gamma,ionQ,f,u,p,s,cv,dpr,dpT
    integer, intent(out) :: phase,err
    real(dp), dimension(ncharged) :: Z, A
    real(dp) :: Mu
    real(dp) :: rsi,ionQ2,fi,ui,si,pi,cvi,dpri,dpti,fm,um,sm,pm,cvm,dprm,dptm
    integer :: i
    
    err = 0
    Mu = amu/Melectron
    Z = real(nuclib% Z(charged_ids(1:ncharged)),dp)
    A = real(nuclib% A(charged_ids(1:ncharged)),dp)
    
    ! set the phase
    rsi = rs*ionic% Z73*Mu
    Gamma = Gamma_e * ionic% Z53
!     if (Gamma < rq% Gamma_melt .or. rsi < rq% rsi_melt) then
    if (Gamma < rq% Gamma_melt) then
        phase = liquid_phase
    else
        phase = solid_phase
    end if

    ! estimate T_p/T for a mixture and check that quantum effects are not important in liquid phase
!   ionQ2 = 3.0*Gamma_e**2/rs/Mu * sum(Zmix**2*Ymix/Amix) *ionic%A / ionic%Z
    ionQ2 = 3.0*Gamma_e**2/rs/Mu * ionic% Z2XoA2 *ionic% A / ionic% Z
    ionQ = sqrt(ionQ2)
    if (phase == liquid_phase .and. ionQ2 > Q2_threshold) then
        err = strong_quantum_effects
        if (.not. rq% suppress_warnings) &
        &   write(error_unit,*) 'ion_mixture: strong quantum effects in liquid encountered'
    end if

    ! sum contributions
    f = 0.0; u = 0.0; p = 0.0; s = 0.0; cv = 0.0; dpr = 0.0; dpT = 0.0
    do i = 1, ncharged
        if (Yion(i) > rq% Ythresh) then
            call one_ion(Gamma_e,rs,Z(i),A(i),phase,fi,ui,pi,si,cvi,dpri,dpti)
            f = f + Yion(i)*fi
            u = u + Yion(i)*ui
            s = s + Yion(i)*(si-log(Yion(i)))
            p = p + Yion(i)*pi
            cv = cv + Yion(i)*cvi
            dpr = dpr + Yion(i)*dpri
            dpt = dpt + Yion(i)*dpti
        end if
    end do
    call LM_corrections(rs,Gamma_e,ionic,fm,um,pm,sm,cvm,dprm,dptm)
    f = f*ionic%A + fm
    u = u*ionic%A + um
    s = s*ionic%A + sm
    p = p*ionic%A + pm
    cv = cv*ionic%A + cvm
    dpr = dpr*ionic%A + dprm
    dpt = dpt*ionic%A + dptm
end subroutine ion_mixture

subroutine LM_corrections(rs,Gamma_e,ionic,f,u,p,s,cv,dpr,dpt)
    use constants_def
    use nucchem_def
    real(dp), intent(in) :: rs,Gamma_e
    type(composition_info_type), intent(in) :: ionic
    real(dp), intent(out) :: f,u,p,s,cv,dpr,dpt
    real(dp), parameter :: threshold = 1.0e-9
    real(dp) :: Gamma, dif0, difR,difFDH,D,p3,d0,gp,f0,q,r,gq,g,g_dg,u_dg
    
    Gamma = Gamma_e*ionic%Z53
    if (rs < threshold) then
        dif0 = ionic%Z52-sqrt(ionic%Z2**3/ionic%Z)
    else
        dif0 = ionic%ZZ1_32 - sqrt((ionic%Z2+ionic%Z)**3/ionic%Z)
    end if
    difR = dif0/ionic%Z52
    difFDH = dif0*Gamma_e*sqrt(onethird*Gamma_e)
    D = ionic%Z2/ionic%Z**2
    
    f = 0.0; u = 0.0; p = 0.0; s = 0.0; cv = 0.0; dpr = 0.0; dpt = 0.0
    if (abs(D-1.0) < threshold) return  ! no correction necessary
    p3 = D**(-0.2)
    d0 = (2.6*difR + 14.0*difR**3)/(1.0-p3)
    gp = d0*Gamma**p3
    f0 = difFDH/(1.0+gp)
    q = D**2 * 0.0117
    r = 1.5/p3 - 1.0
    gq = q*gp
    f = f0/(1.0+gq)**r
    g = 1.5-p3*gp/(1.0+gp)-r*p3*gq/(1.0+gq)
    u = f*g
    s = u-f
    p = onethird*u
    g_dg = -p3**2*(gp/(1.0+gp)**2 + r*gq/(1.0+gq)**2)
    u_dg = u*g + f*g_dg
    cv = u-u_dg
    dpr = p+u_dg/9.0
    dpt = p-onethird*u_dg
end subroutine LM_corrections

subroutine one_ion(Gamma_e,rs,Z,A,phase,f,u,p,s,cv,dpr,dpt)
    use iso_fortran_env, only: error_unit
    use constants_def
    use dStar_eos_def

    real(dp), intent(in) :: Gamma_e,rs,Z,A
    integer, intent(in) :: phase
    real(dp), intent(out) :: f,u,p,s,cv,dpr,dpt
    real(dp) :: Gamma,theta
    real(dp) :: fie,uie,pie,sie,cvie,dprie,dptie
    real(dp) :: fii,uii,pii,sii,cvii,dprii,dptii
    real(dp), dimension(number_lattice_products) :: fiil, uiil, piil, siil, cviil, dpriil, dptiil
    real(dp), parameter :: sevensixth = 7.0/6.0
    Gamma = Gamma_e*Z**fivethird
    theta = Gamma*sqrt(3.0*Melectron/A/amu/rs)/Z**sevensixth
    select case (phase)
        case (liquid_phase)
            call ie_liquid(Gamma_e,rs,Z,fie,uie,pie,sie,cvie,dprie,dptie)
            call ii_liquid(Gamma,theta,fii,uii,pii,sii,cvii,dprii,dptii)
        case (solid_phase)
            call ie_lattice(Gamma,theta,rs,Z,fie,uie,pie,sie,cvie,dprie,dptie)
            call ii_lattice(Gamma,theta,fiil,uiil,piil,siil,cviil,dpriil,dptiil)
            fii = fiil(total); uii = uiil(total); pii = piil(total); sii = siil(total)
            cvii = cviil(total); dprii = dpriil(total); dptii = dptiil(total)
        case default
            fie = 0.0; fii = 0.0; pie = 0.0; pii = 0.0
            sie = 0.0; sii = 0.0; cvie = 0.0; cvii = 0.0
            dprie = 0.0; dprii = 0.0; dptie = 0.0; dptii = 0.0
            write (error_unit,'(a)') 'ion_eos: one_ion:: bad value of ion phase'
    end select
    
    f = fie + fii
    u = uie + uii
    p = pie + pii
    s = sie + sii
    cv = cvie + cvii
    dpr = dprie + dprii
    dpt = dptie + dptii
end subroutine one_ion

subroutine ie_liquid(Gamma_e,rs,Z,f,u,p,s,cv,dpr,dpt)
    use constants_def
    
    real(dp), intent(in) :: Gamma_e, rs, Z
    real(dp), intent(out) :: f,u,p,s,cv,dpr,dpt
    real(dp), parameter :: s3 = sqrt(3.0), CTFfac=(18.0/175.0)*(12.0/pi)**(2.0*onethird)
    real(dp) :: Z13,cDH,cTF,xr,xrs,lZ
    real(dp) :: a,b,nu,sGam,nuGam,nuGam_g,nuGam_gg,ty1,ty1_g,ty1_gg,ty2,ty2_x,ty2_xx
    real(dp) :: ty,tyx,ty_x,ty_g,p1,cor1,cor1_x,cor1_g,cor1_xx,cor1_gg,cor1_xg
    real(dp) :: u0,u0_x,u0_g,u0_xx,u0_gg,u0_xg,d0,d0_g,d0_x,d0_xx
    real(dp) :: cor0,cor0_x,cor0_g,cor0_xx,cor0_gg,cor0_xg
    real(dp) :: rgam,q1,q2,h1u,h1d,h1,h1x,h1_x,h1_xx,up,up_x,up_g,up_xx,up_gg,up_xg
    real(dp) :: dn1,dn1_x,dn1_g,dn1_xx,dn1_gg,dn1_xg,dn,dn_x,dn_xx,dn_g,dn_gg,dn_xg
    real(dp) :: fx,fx_g,f_x,fg,f_g,f_gh,f_xx,f_gg,f_xg
    
    Z13 = Z**onethird
    cDH = (Z/s3)*((1.0+Z)**1.5 - 1.0 - Z**1.5)
    cTF = CTFfac*Z**2*(Z13-1.0 + 0.2/sqrt(Z13))
    xrs = finestructure*(9.0*pi/4.0)**onethird
    xr = xrs/rs
    lZ = log(Z)

    a = 1.11*Z**0.475
    b = 0.2+0.078*lZ**2
    nu = 1.16 + 0.08*lZ
    sGam = sqrt(Gamma_e)
    nuGam = Gamma_e**nu
    nuGam_g = nu*nuGam/Gamma_e
    nuGam_gg = (nu-1.0)*nuGam_g/Gamma_e
    ty1 = 1.0/(0.001*Z**2+2.0*Gamma_e)
    ty1_g = -2.0*ty1**2
    ty1_gg = -4.0*ty1*ty1_g
    ty2 = 1.0+6.0*rs**2
    ty2_x = -12.0*rs**2/xr
    ty2_xx = -3.0*ty2_x/xr
    ty = rs**3/ty2*(1.0+ty1)
    tyx = 3.0/xr + ty2_x/ty2
    ty_x = -ty*tyx
    ty_g = rs**3 * ty1_g/ty2
    p1 = (Z-1.0)/9.0
    cor1 = 1.0+p1*ty
    cor1_x = p1*ty_x
    cor1_g = p1*ty_g
    cor1_xx = p1*(ty*(3.0/xr**2 + (ty2_x/ty2)**2-ty2_xx/ty2)-ty_x*tyx)
    cor1_gg = p1*rs**3*ty1_gg/ty2
    cor1_xg = -p1*ty_g*tyx
    u0 = 0.78*sqrt(Gamma_e/Z)*rs**3
    u0_x = -3.0*u0/xr
    u0_g = 0.5*u0/Gamma_e
    u0_xx = -4.0*u0_x/xr
    u0_gg = -0.5*u0_g/Gamma_e
    u0_xg = -3.0*u0_g/xr
    d0_g = Z**3
    d0 = Gamma_e*d0_g + 21.0*rs**3
    d0_x = -63.0*rs**3/xr
    d0_xx = 252.0*rs**3/xr**2
    cor0 = 1.0 + u0/d0
    cor0_x = (u0_x-u0*d0_x/d0)/d0
    cor0_g = (u0_g-u0*d0_g/d0)/d0
    cor0_xx = (u0_xx-(2.0*u0_x*d0_x+u0*d0_xx)/d0+2.0*(d0_x/d0)**2)/d0
    cor0_gg = (u0_gg-2.0*u0_g*d0_g/d0 + 2.0*u0*(d0_g/d0)**2)/d0
    cor0_xg = (u0_xg-(u0_x*d0_g+u0_g*d0_x)/d0 + 2.0*u0*d0_x*d0_g/d0**2)/d0
    
    rgam = sqrt(1.0+xr**2)
  q1 = 0.18/z**0.25
  q2 = 0.2+0.37/sqrt(z)
  h1u = 1.0 + 0.2*xr**2
  h1d = 1.0 + q1*xr + q2*xr**2
  h1=h1u/h1d
  h1x = 0.4*xr/h1u-(q1+2.*q2*xr)/h1d
  h1_x = h1*h1x
  h1_xx = h1_x*h1x+ h1*(.4/h1u-(.4*xr/h1u)**2-2.*q2/h1d+((q1+2.*q2*xr)/h1d)**2)
  up = cDH*sGam+a*cTF*nuGam*cor0*h1
  up_x = a*cTF*nuGam*(cor0_x*h1+cor0*h1_x)
  up_g = 0.5*cDH/sGam+a*cTF*(nuGam_g*cor0+nuGam*cor0_g)*h1
  up_xx = a*cTF*nuGam*(cor0_xx*h1+2.*cor0_x*h1_x+cor0*h1_xx)
  up_gg = -0.25*cDH/(sGam*Gamma_e)  &
    & + a*cTF*(nuGam_gg*cor0+2.0*nuGam_g*cor0_g+nuGam*cor0_gg)*h1
  up_xg = a*cTF*(nuGam_g*(cor0_x*h1+cor0*h1_x) + nuGam*(cor0_xg*h1+cor0_g*h1_x))
  dn1 = b*sGam+a/rs*nuGam*cor1
  dn1_x = a*nuGam*(cor1/xrs+cor1_x/rs)
  dn1_g = 0.5*b/sGam+a/rs*(nuGam_g*cor1+nuGam*cor1_g)
  dn1_xx = a*nuGam/xrs*(2.*cor1_x+xr*cor1_xx)
  dn1_gg = -0.25*b/(Gamma_e*sGam) + a/rs*(nuGam_gg*cor1+2.*nuGam_g*cor1_g+nuGam*cor1_gg)
  dn1_xg = a*(nuGam_g*(cor1/xrs+cor1_x/rs)+nuGam*(cor1_g/xrs+cor1_xg/rs))
  dn = 1.0 + dn1/rgam
  dn_x = dn1_x/rgam-xr*dn1/rgam**3
  dn_xx = (dn1_xx-((2.*xr*dn1_x+dn1)-3.*xr**2*dn1/rgam**2)/rgam**2)/rgam
  dn_g = dn1_g/rgam
  dn_gg = dn1_gg/rgam
  dn_xg = dn1_xg/rgam-xr*dn1_g/rgam**3
  f = -up/dn*Gamma_e
  fx = (up*dn_x/dn-up_x)/dn
  fx_g = ((up_g*dn_x+up_x*dn_g+up*dn_xg-2.0*up*dn_x*dn_g/dn)/dn - up_xg)/dn
  f_x = fx*Gamma_e
  fg = (up*dn_g/dn-up_g)/dn
  f_g = fg*Gamma_e-up/dn
  f_xx = ((up*dn_xx+2.0*(up_x*dn_x-up*dn_x**2/dn))/dn-up_xx)/dn*Gamma_e
  f_gg = 2.0*fg+Gamma_e*((2.*dn_g*(up_g-up*dn_g/dn)+up*dn_gg)/dn-up_gg)/dn
  f_xg = fx+Gamma_e*fx_g
  u = Gamma_e*f_g
    s = u-f
  cv = -Gamma_e**2*f_gg
  p = onethird*(xr*f_x+Gamma_e*f_g)
  dpt = -Gamma_e**2*(xr*fx_g+f_gg)*onethird
  dpr = (12.0*p+xr**2*f_xx+2.0*xr*Gamma_e*f_xg+Gamma_e**2*f_gg)/9.0
end subroutine ie_liquid

subroutine ie_lattice(Gamma,theta,rs,Z,f,u,p,s,cv,dpr,dpt)
    use constants_def
    
    real(dp), intent(in) :: Gamma,theta,rs,Z
    real(dp), intent(out) :: f,u,p,s,cv,dpr,dpt
    real(dp), parameter, dimension(4) :: a = [1.1866,0.684,17.9,41.5]
    real(dp), parameter :: e1 = exp(1.0), px = 0.205
    real(dp) :: aTF
    real(dp) :: xrs,xr,lZ,r3,b(4),finf,finfx,finf_x,finf_xx,q1u,q1d,q1,q1x,q1x_x,q1_x,q1_xx
    real(dp) :: y0,y0_x,y0_g,y0_xx,y0_gg,y0_xg,y1,y1_x,y1_g,y1_xx,y1_gg,y1_xg
    real(dp) :: sa,supa,supa_x,supa_g,supa_xx,supa_gg,supa_xg
    real(dp) :: em2,sb,em2y1,supb,supb_x,supb_g,supb_xx,supb_gg,supb_xg
    real(dp) :: sup,supx,sup_x,supg,sup_g,sup_xx,sup_gg,sup_xg
    real(dp) :: gr3,gr3x,gr3_x,gr3_xx,gr3g,gr3_g,gr3_gg,gr3_xg,w,w_x,w_g,w_xx,w_gg,w_xg
    real(dp) :: f_x,f_g,f_xx,f_gg,f_xg

    aTF = (54.0/175.0)*(12.0/pi)**onethird*finestructure
    xrs = finestructure*(9.0*pi/4.0)**onethird
    xr = xrs/rs

    lZ = log(Z)
    r3 = 1.0/(1.0+0.01*lZ**1.5+0.097/Z**2)
    b(1) = 1.0-a(1)/Z**0.267 + 0.27/Z
    b(2) = 1.0 + 2.25/Z**onethird*(1.0+a(2)*Z**5+0.222*Z**6)/(1.0+0.222*Z**6)
    b(3) = a(4)/(1.0+lZ)
    b(4) = 0.395*lZ + 0.347/Z**1.5
    
    finf = aTF*Z**twothird*b(1)*sqrt(1.0+b(2)/xr**2)
    finfx = -b(2)/((b(2)+xr**2)*xr)
    finf_x = finf*finfx
    finf_xx = finf_x*finfx - finf_x*(b(2)+3.0*xr**2)/((b(2)+xr**2)*xr)
    
    q1u = b(3)+a(3)*xr**2
    q1d = 1.0 + b(4)*xr**2
    q1 = q1u/q1d
    q1x = 2.0*xr*(a(3)/q1u - b(4)/q1d)
    q1x_x = q1x/xr + 4.0*xr**2*((b(4)/q1d)**2 - (a(3)/q1u)**2)
    q1_x = q1*q1x
    q1_xx = q1_x*q1x + q1*q1x_x

    if (theta < 6.0/px) then
        y0 = (px*theta)**2
        y0_x = y0/xr
        y0_g = 2.0*y0/Gamma
        y0_xx = 0.0
        y0_gg = y0_g/Gamma
        y0_xg = y0_g/xr
        y1 = exp(y0)
        y1_x = y1*y0_x
        y1_g = y1*y0_g
        y1_xx = y1*(y0_x**2+y0_xx)
        y1_gg = y1*(y0_g**2+y0_gg)
        y1_xg = y1*(y0_x*y0_g+y0_xg)
        sa = 1.0+y1
        supa = log(sa)
        supa_x = y1_x/sa
        supa_g = y1_g/sa
        supa_xx = (y1_xx-y1_x**2/sa)/sa
        supa_gg = (y1_gg-y1_g**2/sa)/sa
        supa_xg = (y1_xg-y1_x*y1_g/sa)/sa
        em2 = e1-2.0
        sb = e1-em2/y1
        supb = log(sb)
        em2y1 = em2/(y1**2*sb)
        supb_x = em2y1*y1_x
        supb_g = em2y1*y1_g
        supb_xx = em2y1*(y1_xx-2.0*y1_x**2/y1-y1_x*supb_x)
        supb_gg = em2y1*(y1_gg-2.0*y1_g**2/y1-y1_g*supb_g)
        supb_xg = em2y1*(y1_xg-2.0*y1_x*y1_g/y1-y1_g*supb_x)
        sup = sqrt(supa/supb)
        supx = 0.5*(supa_x/supa-supb_x/supb)
        sup_x = sup*supx
        supg = 0.5*(supa_g/supa-supb_g/supb)
        sup_g = sup*supg
        sup_xx = sup_x*supx + sup*0.5*(supa_xx/supa-(supa_x/supa)**2 - supb_xx/supb+(supb_x/supb)**2)
        sup_gg = sup_g*supg + sup*0.5*(supa_gg/supa-(supa_g/supa)**2 - supb_gg/supb+(supb_g/supb)**2)
        sup_xg = sup_x*supg + sup*0.5*((supa_xg-supa_x*supa_g/supa)/supa  &
                & - (supb_xg-supb_x*supb_g/supb)/supb)
    else
        sup=px*theta
        sup_x = 0.5*px*theta/xr
        sup_g = px*theta/Gamma
        sup_xx = -0.5*sup_x/xr
        sup_gg = 0.0
        sup_xg = sup_x/Gamma
    end if
    gr3 = (Gamma/sup)**r3
    gr3x = -r3*sup_x/sup
    gr3_x = gr3*gr3x
    gr3_xx = gr3_x*gr3x - r3*gr3*(sup_xx/sup-(sup_x/sup)**2)
    gr3g = r3*(1.0/Gamma-sup_g/sup)
    gr3_g = gr3*gr3g
    gr3_gg = gr3_g*gr3g + gr3*r3*((sup_g/sup)**2-sup_gg/sup-1.0/Gamma**2)
    gr3_xg = gr3_g*gr3x+gr3*r3*(sup_x*sup_g/sup**2 - sup_xg/sup)
    w = 1.0+q1/gr3
    w_x = q1_x/gr3-q1*gr3_x/gr3**2
    w_g = -q1*gr3_g/gr3**2
    w_xx = q1_xx/gr3-(2.0*q1_x*gr3_x+q1*(gr3_xx-2.0*gr3_x**2/gr3))/gr3**2
    w_gg = q1*(2.0*gr3_g**2/gr3-gr3_gg)/gr3**2
    w_xg = -(q1_x*gr3_g+q1*(gr3_xg-2.0*gr3_x*gr3_g/gr3))/gr3**2
    
    f = -Gamma*finf*w
    f_x = -Gamma*(finf_x*w+finf*w_x)
    f_xx = -Gamma*(finf_xx*w + 2.0*finf_x*w_x + finf*w_xx)
    f_g = -finf*w-Gamma*finf*w_g
    f_gg = -2.0*finf*w_g - Gamma*finf*w_gg
    f_xg = -finf_x*w-finf*w_x-Gamma*(finf_x*w_g+finf*w_xg)
    
    s = -Gamma**2*finf*w_g
    u = s+f
    cv = -Gamma**2*f_gg
    p = onethird*(xr*f_x+Gamma*f_g)
    dpt = onethird*Gamma**2*(xr*finf*(finfx*w_g+w_xg)-f_gg)
    dpr = (12.0*p+xr**2*f_xx + 2.0*xr*Gamma*f_xg + Gamma**2*f_gg)/9.0
end subroutine ie_lattice

subroutine ii_liquid(Gamma,theta,f,u,p,s,cv,dpr,dpt)
    use constants_def
    
    real(dp), intent(in) :: Gamma, theta
    real(dp), intent(out) :: f, u, p, s, cv, dpr, dpt
    real(dp), parameter, dimension(3) :: A = [-0.907347,0.62849, &
            & -0.5*sqrt(3.0)+ 0.907347/sqrt(0.62849)]
    real(dp), parameter, dimension(4) :: B = [4.50e-3, 170.0, -8.4e-5, 3.70e-3]
    real(dp) :: fq, f0, fid
    
    f0 = A(1)*(sqrt(Gamma*(A(2)+Gamma))-A(2)*log(sqrt(Gamma/A(2)) &
            & +sqrt(1.0+Gamma/A(2))))  &
            & + 2.0*A(3)*(sqrt(Gamma)-atan(sqrt(Gamma)))
    
    f = f0 + B(1)*(Gamma-B(2)*log(1.0+Gamma/B(2))) + 0.5*B(3)*log(1.0+Gamma**2/B(4))
    u = Gamma**1.5*(A(1)/sqrt(A(2)+Gamma) + A(3)/(1.0+Gamma)) +  &
            & Gamma**2 * (B(1)/(B(2)+Gamma) + B(3)/(B(4)+Gamma**2))
    cv = 0.5*Gamma**1.5*(A(3)*(Gamma-1.0)/(Gamma+1.0)**2  &
            & - A(1)*A(2)/(Gamma+A(2))**1.5)  &
            & + Gamma**2*(B(3)*(Gamma**2-B(4))/(Gamma**2+B(4))**2  &
            & - B(1)*B(2)/(Gamma+B(2))**2)
    p = onethird*u
    dpr = (4.0*u-cv)/9.0
    dpt = onethird*cv
    
    ! ideal component
    fid = 3.0*log(theta) - 1.5*log(Gamma) - 0.5*log(6.0/pi) - 1.0
    f = f + fid
    u = u + 1.5
    p = p + 1.0
    cv = cv + 1.5
    s = u-f
    dpt = dpt + 1.0
    dpr = dpr + 1.0
    
    ! addition of thermal component, Potekhin & Chabrier (2010)
    fq = theta**2/24.0
    f = f + fq
    u = u + 2.0*fq
    s = s + fq
    cv = cv - 2.0*fq
    p = p + fq
    dpr = dpr + 2.0*fq
    dpt = dpt - fq
end subroutine ii_liquid

subroutine ii_lattice(Gamma,theta,f,u,p,s,cv,dpr,dpt)
    use constants_def
    
    real(dp), intent(in) :: Gamma, theta
    real(dp), intent(out), dimension(number_lattice_products) :: f, u, p, s, cv, dpr, dpt
    real(dp), parameter :: KM = -0.895929255682, w=0.5113875, clm = -2.49389
    real(dp), parameter, dimension(0:8) :: alpha = &
                    & [0.0,0.932446,0.334547,0.265764,0.0,0.0,4.757014e-3,0.0,4.7770935e-3]
    real(dp), parameter, dimension(0:8) :: a =  &
                    & [1.0,0.1839,0.593586,5.4814e-3,5.01813e-4,0.0,3.9247e-7,0.0,5.8356e-11]
    real(dp), parameter, dimension(0:7) :: b =  &
                    & [261.66,0.0,7.07997,0.0,0.0409484,3.97355e-4,5.11148e-5,2.19749e-6]
    real(dp), parameter, dimension(3) :: afh = [10.9,247.0,1.765e5]
    real(dp), parameter :: b1 = 0.12, c1 = b1/afh(1)
    real(dp), parameter :: theta_low = 1.0e-5, theta_high = 1.0e5
    real(dp), dimension(0:2) :: Ac,Bc
    real(dp) :: ea, omea, e0, madelung, c1t2, sup, tq, theta2, theta4, gninv, ninv, pfac
    integer :: n
        
    f = 0.0; u = 0.0; p = 0.0; s = 0.0; cv = 0.0; dpt = 0.0; dpr = 0.0

    ! harmonic component
    if (theta > theta_high) then
        u(harmonic) = 3.0/(alpha(6)*theta**3)
        f(harmonic) = -u(harmonic)/3.0
        cv(harmonic) = 4.0*u(harmonic)
    else if (theta < theta_low) then
        f(harmonic) = 3.0*log(theta) + clm - 1.5*w*theta + theta**2/24.0
        u(harmonic) = 3.0 - 1.5*w*theta + theta**2/12.0
        cv(harmonic) = 3.0 - theta**2/12.0
    else
        Ac = Ath(theta)
        Bc = Bth(theta)
        do n = 1,3
            ea = exp(-alpha(n)*theta)
            omea = 1.0-ea
            f(harmonic) = f(harmonic) + log(omea)
            u(harmonic) = u(harmonic) + alpha(n)*ea/omea
            cv(harmonic) = cv(harmonic) + alpha(n)**2*ea/omea**2
        end do
        f(harmonic) = f(harmonic) - Ac(0)/Bc(0)
        u(harmonic) = (u(harmonic)-(Ac(1)*Bc(0)-Ac(0)*Bc(1))/Bc(0)**2)*theta
        cv(harmonic) = cv(harmonic) + &
        & ((Ac(2)*Bc(0) - Bc(2)*Ac(0))*Bc(0)  &
        & - 2.0*(Ac(1)*Bc(0)-Bc(1)*Ac(0))*Bc(1))/Bc(0)**3
        cv(harmonic) = cv(harmonic)*theta**2
    end if
    e0 = 1.5*w*theta
    s(harmonic) = u(harmonic) - f(harmonic)
    u(harmonic) = u(harmonic) + e0
    f(harmonic) = f(harmonic) + e0
    p(harmonic) = 0.5*u(harmonic)
    dpt(harmonic) = 0.5*cv(harmonic)
    dpr(harmonic) = 0.75*u(harmonic) - 0.25*cv(harmonic)
    
    ! anharmonic component
    theta2 = theta**2
    theta4 = theta**4
    c1t2 = c1*theta2
    sup = exp(-c1t2)
    tq = b1*theta2/Gamma
    gninv = sup
    do n = 1,3
        gninv = gninv/Gamma
        ninv = 1.0/n
        f(anharmonic) = f(anharmonic) - ninv*afh(n)*gninv
        u(anharmonic) = u(anharmonic) + gninv*afh(n)*(1.0+2.0*ninv*c1t2)
        pfac = onethird*afh(n) + c1t2*afh(n)*ninv
        p(anharmonic) = p(anharmonic) + pfac*gninv
        cv(anharmonic) = cv(anharmonic) + gninv*((1.0+n)*afh(n)  &
                & + (4.0-2.0*ninv)*afh(n)*c1t2 + 4.0*afh(n)*c1**2*ninv*theta4)
        dpt(anharmonic) = dpt(anharmonic) + (pfac*(1.0+n+2.0*c1t2)-2.0*ninv*afh(n)*c1t2)*gninv
        dpr(anharmonic) = dpr(anharmonic) + (pfac*(1.0-onethird*n-c1t2)+afh(n)*ninv*c1t2)*gninv
    end do
    f(anharmonic) = f(anharmonic) - tq
    u(anharmonic) = u(anharmonic) - tq
    p(anharmonic) = p(anharmonic) - 2.0*onethird*tq
    dpr(anharmonic) = dpr(anharmonic) - tq/4.5
    
    ! now add the Madelung and zero-point energies
    madelung = KM*Gamma
    f(total) = f(harmonic) + f(anharmonic) + madelung
    u(total) = u(harmonic) + u(anharmonic) + madelung
    p(total) = p(harmonic) + p(anharmonic) + onethird*madelung
    s(total) = s(harmonic)
    cv(total) = cv(harmonic) + cv(anharmonic)
    dpt(total) = dpt(harmonic) + dpt(anharmonic)
    dpr(total) = dpr(harmonic) + dpr(anharmonic) + madelung/2.25
    
    contains
    function Ath(theta)
        real(dp), intent(in) :: theta
        real(dp), dimension(0:2) :: Ath
        integer :: n
        Ath(0) = a(8)
        Ath(1) = 0.0
        Ath(2) = 0.0
        do n = 7,0,-1
            Ath(2) = Ath(2)*theta+Ath(1)
            Ath(1) = Ath(1)*theta+Ath(0)
            Ath(0) = Ath(0)*theta + a(n)
        end do
        Ath(2) = Ath(2) * 2.0
    end function Ath
    function Bth(theta)
        real(dp), intent(in) :: theta
        real(dp), dimension(0:2) :: Bth
        integer :: n
        Bth(0) = b(7)
        Bth(1) = 0.0
        Bth(2) = 0.0
        do n = 6,0,-1
            Bth(2) = Bth(2)*theta+Bth(1)
            Bth(1) = Bth(1)*theta+Bth(0)
            Bth(0) = Bth(0)*theta + b(n)
        end do
        Bth(2) = Bth(2)*2.0
        Bth(0) = Bth(0) + alpha(6)*a(6)*theta**9 + alpha(8)*a(8)*theta**11
        Bth(1) = Bth(1) + 9.0*alpha(6)*a(6)*theta**8 + 11.0*alpha(8)*a(8)*theta**10
        Bth(2) = Bth(2) + 72.0*alpha(6)*a(6)*theta**7 + 110.0*alpha(8)*a(8)*theta**9
    end function Bth
end subroutine ii_lattice
end module ion_eos
