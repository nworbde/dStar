! This code is a customized version of the routine for electron exchange interactions made availiby by
! A.Y. Potekhin at http://www.ioffe.ru/astro/EIP
!
!

module electron_eos

contains

subroutine get_helm_eos_results(rho,T,ionic,f,u,p,s,cv,dpr,dpt,eta)
    use constants_def, only: dp
    use dStar_eos_def
    use helm
    use nucchem_def, only: composition_info_type
    real(dp), intent(in) :: rho,T
    type(composition_info_type), intent(in) :: ionic
    real(dp), intent(out) :: f,u,p,s,cv,dpr,dpt,eta
    real(dp) :: logT,logRho,Zfrac,Xfrac,abar,zbar
    real(dp), dimension(num_helm_results) :: helm_res
    logical, parameter :: clip = .true.
    integer :: ierr
    
    ! Xfrac and Zfrac are not used by helm, so set to a nonsensical value
    Xfrac = -1.0; Zfrac = -1.0
    f = 0.0; u = 0.0; p = 0.0; s = 0.0; cv = 0.0; dpr = 0.0; dpt = 0.0

    ! check for pure neutron matter
    if (abs(ionic% Yn - 1.0) <= epsilon(1.0)) return
    
    ! set inputs and call helm
    abar = ionic% A/(1.0-ionic% Yn)
    zbar = ionic% Z
    logT = log10(T)
    logRho = log10(rho)
    call helmeos2(T,logT,rho,logRho,Zfrac,Xfrac,abar,zbar,helm_res,clip,ierr)
    if (ierr /= 0) return
    
    ! store the results
    u = helm_res(h_eele)
    p = helm_res(h_pele)
    s = helm_res(h_sele)
    f = u-T*s
    cv = helm_res(h_deept)
    dpr = rho*helm_res(h_dpepd)
    dpt = T*helm_res(h_dpept)
    eta = helm_res(h_etaele)
end subroutine get_helm_eos_results

subroutine ee_exchange(Gamma_e,rs,f,u,p,s,cv,dpr,dpt)
    use constants_def
    real(dp), intent(in) :: Gamma_e,rs
    real(dp), intent(out) :: f,u,p,s,cv,dpr,dpt
    real(dp), parameter :: theta_fac = 2.0*(4.0/9.0/pi)**(2.0/3.0), theta_low = 0.005
    real(dp) :: theta,sqth,theta2,theta3,theta4,cht1,sht1,cht2,sht2
    real(dp) :: t1,t1_h,t1_hh,t2,t2_h,t2_hh,a0,a0_h,a0_hh,a1,a1_h,a1_hh,a,ah,a_h,a_hh
    real(dp) :: b0,b0_h,b0_hh,b1,b1_h,b1_hh,b,bh,b_h,b_hh,d0,d0_h,d0_hh,d1,d1_h,d1_hh,d,dh,d_h,d_hh
    real(dp) :: e0,e0_h,e0_hh,e1,e1_h,e1_hh,e,eh,e_h,e_hh,exp1th,c,c_h,c_hh,discr,di_h,di_hh
    real(dp) :: s1,s1h,s1_h,s1_hh,s1_g,s1_hg,b2,b2_h,b2_hh
    real(dp) :: sqge,s2,s2h,s2_h,s2_hh,s2_g,s2_gg,s2_hg,r3,r3_h,r3_hh,r3_g,r3_gg,r3_hg
    real(dp) :: b3,b3_h,b3_hh,c3,c3_h,c3_hh,s3,s3_h,s3_hh,s3_g,s3_gg,s3_hg,b4,b4_h,b4_hh
    real(dp) :: c4,c4_h,c4_hh,c4_g,c4_gg,c4_hg,s4a,s4ah,s4a_h,s4a_hh,s4b,s4b_h,s4b_hh
    real(dp) :: s4c,s4c_h,s4c_hh,s4c_g,s4c_gg,s4c_hg,s4,s4_h,s4_hh,s4_g,s4_gg,s4_hg
    real(dp) :: up1,dn1,up2,dn2
    real(dp) :: f_h,f_g,f_hh,f_gg,f_hg,pdlh,pdlg
    
    theta = theta_fac*rs/Gamma_e ! non-relativistic degeneracy parameter
    sqth = sqrt(theta)
    theta2 = theta**2
    theta3 = theta2*theta
    theta4 = theta3*theta
    if (theta > theta_low) then
        cht1 = cosh(1.0/theta)
        sht1 = sinh(1.0/theta)
        cht2 = cosh(1.0/sqth)
        sht2 = sinh(1.0/sqth)
        t1 = sht1/cht1 ! tanh(1.d0/theta)
        t2 = sht2/cht2 ! tanh(1./dsqrt(theta))
        t1_h = -1.0/(theta*cht1)**2 ! d t1 / d\theta
        t1_hh = 2.0/(theta*cht1)**3*(cht1-sht1/theta)
        t2_h = -0.5*sqth/(theta*cht2)**2
        t2_hh = (0.75*sqth*cht2-0.5*sht2)/(theta*cht2)**3
    else
        t1 = 1.0
        t2 = 1.0
        t1_h = 0.0
        t2_h = 0.0
        t1_hh = 0.0
        t2_hh = 0.0
    end if
    a0 = 0.75+3.04363*theta2-0.09227*theta3+1.7035*theta4
    a0_h = 6.08726*theta-0.27681*theta2+6.814*theta3
    a0_hh = 6.08726-0.55362*theta+20.442*theta2
    a1 = 1.+8.31051*theta2+5.1105*theta4
    a1_h = 16.62102*theta+20.442*theta3
    a1_hh = 16.62102+61.326*theta2
    a = 0.610887*a0/a1*t1 ! hf fit of perrot and dharma-wardana
    ah = a0_h/a0-a1_h/a1+t1_h/t1
    a_h = a*ah
    a_hh=a_h*ah+a*(a0_hh/a0-(a0_h/a0)**2-a1_hh/a1+(a1_h/a1)**2+t1_hh/t1-(t1_h/t1)**2)
    
    b0 = 0.341308+12.070873*theta2+1.148889*theta4
    b0_h = 24.141746*theta+4.595556*theta3
    b0_hh = 24.141746+13.786668*theta2
    b1 = 1.0+10.495346*theta2+1.326623*theta4
    b1_h = 20.990692*theta+5.306492*theta3
    b1_hh = 20.990692+15.919476*theta2
    b = sqth*t2*b0/b1
    bh = 0.5/theta+t2_h/t2+b0_h/b0-b1_h/b1
    b_h = b*bh
    b_hh = b_h*bh+b*(-0.5/theta2+t2_hh/t2-(t2_h/t2)**2  &
                & +  b0_hh/b0 - (b0_h/b0)**2-b1_hh/b1+(b1_h/b1)**2)

    d0 = 0.614925+16.996055*theta2+1.489056*theta4
    d0_h = 33.99211*theta+5.956224*theta3
    d0_hh = 33.99211+17.868672*theta2
    d1 = 1.0+10.10935*theta2+1.22184*theta4
    d1_h = 20.2187*theta+4.88736*theta3
    d1_hh = 20.2187+14.66208*theta2
    d = sqth*t2*d0/d1
    dh = 0.5/theta+t2_h/t2+d0_h/d0-d1_h/d1
    d_h = d*dh
    d_hh=d_h*dh+d*(-0.5/theta2+t2_hh/t2-(t2_h/t2)**2 &
        & + d0_hh/d0-(d0_h/d0)**2-d1_hh/d1+(d1_h/d1)**2)

    e0 = 0.539409 + 2.522206*theta2+0.178484*theta4
    e0_h = 5.044412*theta+0.713936*theta3
    e0_hh = 5.044412+2.141808*theta2
    e1 = 1.0+2.555501*theta2+0.146319*theta4
    e1_h = 5.111002*theta+0.585276*theta3
    e1_hh = 5.111002+1.755828*theta2
    e = theta*t1*e0/e1
    eh = 1.0/theta+t1_h/t1+e0_h/e0-e1_h/e1
    e_h = e*eh
    e_hh = e_h*eh+e*(t1_hh/t1-(t1_h/t1)**2+e0_hh/e0-(e0_h/e0)**2 &
                & -e1_hh/e1+(e1_h/e1)**2-1.0/theta2)

    exp1th = exp(-1.0/theta)

  c = (0.872496+0.025248*exp1th)*e
    c_h = 0.025248*exp1th/theta2*e+c*e_h/e
    c_hh = 0.025248*exp1th/theta2*(e_h+(1.0-2.0*theta)/theta2*e) &
            & +c_h*e_h/e+c*e_hh/e-c*(e_h/e)**2
    discr = sqrt(4.0*e-d**2)
    di_h = 0.5/discr*(4.0*e_h-2.0*d*d_h)
    di_hh=(-((2.0*e_h-d*d_h)/discr)**2 + 2.0*e_hh-d_h**2-d*d_hh)/discr
    
    s1 = -c/e*Gamma_e
    s1h = c_h/c-e_h/e
    s1_h = s1*s1h
    s1_hh = s1_h*s1h+s1*(c_hh/c-(c_h/c)**2-e_hh/e+(e_h/e)**2)
    s1_g = -c/e ! => s1dgg=0
    s1_hg = s1_g*(c_h/c-e_h/e)
    b2 = b-c*d/e
    b2_h = b_h-(c_h*d+c*d_h)/e+c*d*e_h/e**2
    b2_hh = b_hh-(c_hh*d+2.0*c_h*d_h+c*d_hh)/e  &
            & + (2.0*(c_h*d+c*d_h-c*d*e_h/e)*e_h+c*d*e_hh)/e**2

    sqge = sqrt(Gamma_e)
    s2 = -2.0/e*b2*sqge
    s2h = b2_h/b2-e_h/e
    s2_h = s2*s2h
    s2_hh = s2_h*s2h+s2*(b2_hh/b2-(b2_h/b2)**2-e_hh/e+(e_h/e)**2)
    s2_g = 0.5*s2/Gamma_e
    s2_gg = -0.5*s2_g/Gamma_e
    s2_hg = 0.5*s2_h/Gamma_e

    r3 = e*Gamma_e+d*sqge+1.0
    r3_h = e_h*Gamma_e+d_h*sqge
    r3_hh = e_hh*Gamma_e+d_hh*sqge
    r3_g = e+0.5*d/sqge
    r3_gg = -0.25*d/(Gamma_e*sqge)
    r3_hg = e_h+0.5*d_h/sqge

    b3 = a-c/e
    b3_h = a_h-c_h/e+c*e_h/e**2
    b3_hh = a_hh-c_hh/e+(2.0*c_h*e_h+c*e_hh)/e**2-2.0*c*e_h**2/e**3

    c3 = (d/e*b2-b3)/e ! =d*b2/e**2-b3/e
    c3_h = (d_h*b2+d*b2_h+b3*e_h)/e**2-2.0*d*b2*e_h/e**3-b3_h/e
    c3_hh = (-b3_hh+(d_hh*b2+2.0*d_h*b2_h+d*b2_hh+b3_h*e_h+b3*e_hh+b3_h*e_h)/e &
            & - 2.0*((d_h*b2+d*b2_h+b3*e_h+d_h*b2+d*b2_h)*e_h+d*b2*e_hh)/e**2 &
            & + 6.0*d*b2*e_h**2/e**3)/e

    s3 = c3*log(r3)
    s3_h = s3*c3_h/c3+c3*r3_h/r3
    s3_hh = (s3_h*c3_h+s3*c3_hh)/c3-s3*(c3_h/c3)**2+(c3_h*r3_h+c3*r3_hh)/r3-c3*(r3_h/r3)**2
    s3_g = c3*r3_g/r3
    s3_gg = c3*(r3_gg/r3-(r3_g/r3)**2)
    s3_hg = (c3_h*r3_g+c3*r3_hg)/r3-c3*r3_g*r3_h/r3**2

    b4 = 2.0-d**2/e
    b4_h = e_h*(d/e)**2-2.0*d*d_h/e
    b4_hh = e_hh*(d/e)**2+2.0*e_h*(d/e)**2*(d_h/d-e_h/e)-2.0*(d_h**2+d*d_hh)/e+2.*d*d_h*e_h/e**2

    c4 = 2.0*e*sqge+d
    c4_h = 2.0*e_h*sqge+d_h
    c4_hh = 2.0*e_hh*sqge+d_hh
    c4_g = e/sqge
    c4_gg = -0.5*e/(Gamma_e*sqge)
    c4_hg = e_h/sqge

    s4a = 2.0/e/discr
    s4ah = e_h/e+di_h/discr
    s4a_h = -s4a*s4ah
    s4a_hh = -s4a_h*s4ah-s4a*(e_hh/e-(e_h/e)**2+di_hh/discr-(di_h/discr)**2)
    s4b = d*b3+b4*b2
    s4b_h = d_h*b3+d*b3_h+b4_h*b2+b4*b2_h
    s4b_hh = d_hh*b3+2.0*d_h*b3_h+d*b3_hh+b4_hh*b2+2.*b4_h*b2_h+b4*b2_hh
    s4c = atan(c4/discr)-atan(d/discr)
    up1 = c4_h*discr-c4*di_h
    dn1 = discr**2+c4**2
    up2 = d_h*discr-d*di_h
    dn2 = discr**2+d**2
    s4c_h = up1/dn1-up2/dn2
    s4c_hh = (c4_hh*discr-c4*di_hh)/dn1-up1*2.0*(discr*di_h+c4*c4_h)/dn1**2 &
            & -(d_hh*discr-d*di_hh)/dn2+up2*2.0*(discr*di_h+d*d_h)/dn2**2
    s4c_g = c4_g*discr/dn1
    s4c_gg = c4_gg*discr/dn1-2.0*c4*discr*(c4_g/dn1)**2
    s4c_hg = (c4_hg*discr+c4_g*di_h-c4_g*discr/dn1*2.*(discr*di_h+c4*c4_h))/dn1
    s4 = s4a*s4b*s4c
    s4_h = s4a_h*s4b*s4c+s4a*s4b_h*s4c+s4a*s4b*s4c_h
    s4_hh = s4a_hh*s4b*s4c+s4a*s4b_hh*s4c+s4a*s4b*s4c_hh &
            & + 2.0*(s4a_h*s4b_h*s4c+s4a_h*s4b*s4c_h+s4a*s4b_h*s4c_h)
    s4_g = s4a*s4b*s4c_g
    s4_gg = s4a*s4b*s4c_gg
    s4_hg = s4a*s4b*s4c_hg+s4c_g*(s4a_h*s4b+s4a*s4b_h)

    f = s1+s2+s3+s4
    f_h = s1_h+s2_h+s3_h+s4_h
    f_g = s1_g+s2_g+s3_g+s4_g
    f_hh = s1_hh+s2_hh+s3_hh+s4_hh
    f_gg = s2_gg+s3_gg+s4_gg
    f_hg = s1_hg+s2_hg+s3_hg+s4_hg
    
    p = onethird*(Gamma_e*f_g-2.0*theta*f_h)
    u = Gamma_e*f_g-theta*f_h
    s = (Gamma_e*s2_g-s2+Gamma_e*s3_g-s3+s4a*s4b*(Gamma_e*s4c_g-s4c))-theta*f_h
    if (abs(s) < 1.0e-9*abs(theta*f_h)) s = 0.0 ! accuracy loss
    cv = 2.0*theta*(Gamma_e*f_hg-f_h)-theta**2*f_hh-Gamma_e**2*f_gg
    if (abs(cv) < 1.0e-9*abs(Gamma_e**2*f_gg)) cv = 0.0 ! accuracy
    pdlh = onethird*theta*(Gamma_e*f_hg-2.0*f_h-2.0*theta*f_hh)
    pdlg = onethird*Gamma_e*(f_g+Gamma_e*f_gg-2.0*theta*f_hg)
    dpr = p+onethird*(pdlg-2.0*pdlh)
    dpt = Gamma_e*(theta*f_hg - onethird*Gamma_e*f_gg)-theta*(f_h/0.75+theta*f_hh/1.5)

end subroutine ee_exchange

end module electron_eos
