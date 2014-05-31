module eval_crust_neutrino
    use neutrino_def
    
    ! number of neutrino flavors, in addition to electron neutrino
    real(dp), parameter :: cn_dn = 2.0
    real(dp), parameter :: cn_snthtw = 0.2319
    real(dp), parameter :: cn_cv = 0.5+2.0*cn_snthtw
    real(dp), parameter :: cn_ca = 0.5
    real(dp), parameter :: cn_cvp = 1.0-cn_cv
    real(dp), parameter :: cn_cap = 1.0-cn_ca
    
    contains
    
    subroutine emissivity(rho,T,ionic,chi,Tcn,include_channel,eps)
        use constants_def
        use nucchem_def, only: composition_info_type
        real(dp), intent(in) :: rho,T,chi
        type(composition_info_type), intent(in) :: ionic
        real(dp), intent(in) :: Tcn
        logical, dimension(num_crust_nu_channels) :: include_channel
        type(crust_neutrino_emissivity_channels), intent(out) :: eps
        real(dp) :: Tme
        real(dp) :: rY,lambda, nn, kn       
        
        Tme  = Melectron*clight2/boltzmann
        rY = rho*ionic%Ye
        lambda = T/Tme
        
        call clear_epsilon
        if (include_channel(icrust_nu_pair))  &
                & call pair(T,rY,eps% pair)
        if (include_channel(icrust_nu_photo)) &
                &   call photo(eps% photo,rY,lambda)
        if (include_channel(icrust_nu_plasma)) &
                &   call plasma(T,rY,eps% plasma)
        if (include_channel(icrust_nu_brems)) &
                &   call bremsstrahlung(rho,T,ionic,eps% brems)
        if (include_channel(icrust_nu_pbf)) then
            nn = rho*ionic%Yn/amu/(1.0-chi)
            kn = (1.5*pi**2*nn)**onethird / cm_to_fm
            call pbf(kn,T,Tcn,eps% pbf)
        end if
        eps% total = eps% pair + eps% photo + eps% plasma + eps% brems + eps% pbf
        
        contains
        subroutine clear_epsilon()
            eps% total = 0.0
            eps% pair = 0.0
            eps% photo = 0.0
            eps% plasma = 0.0
            eps% brems = 0.0
            eps% pbf = 0.0
        end subroutine clear_epsilon
    end subroutine emissivity
    
    subroutine plasma(t,r,q)
        ! plasma.f
        ! modified from formulae routine
        ! Itoh et al. ApJSS 102: 411 (1996)
        !
        ! Input:        t = temperature [K]
        !                       r = density*Ye [g cm**-3]
        ! Returns:  q = plasma neutrino emissivity [ergs/g/s]
        !
        real(dp), intent(in) :: t, r
        real(dp), intent(out) :: q
        real(dp) :: rdouble,x,y,enum,fxy,lamb,gamd,gam,ft,fl,qap
        ! equation (4.1)
        real(dp), parameter :: cvd = cn_cv**2 + cn_dn*cn_cvp**2

        ! equations (4.9)-(4.10)
        rdouble = log10(2.0*r)
        x = (  17.5 + rdouble - 3.0*t )/6.0
      y = ( -24.5 + rdouble + 3.0*t )/6.0
        enum = min(0.0,y-1.6+1.25*x)
        if (abs(x)>0.7 .or. y < 0.0) then
            fxy = 1.0
        else
            fxy = 1.05 + (0.39 -1.25*x - 0.35*sin(4.5*x) -0.3*exp(-(4.5*x + 0.9)**2)) &
                                & * exp(-(enum/(0.57-0.25*x))**2)
        end if

        ! equation (2.10)
        lamb = 1.686e-10*t

        ! equations (4.6)-(4.8)
        gamd = 1.1095e11*r/((t**2)*sqrt(1.0+(1.019e-6*r)**(2./3.)))
        gam = sqrt(gamd)
      ft   = 2.4 + 0.6*sqrt(gam) + 0.51*gam + 1.25*gam**1.5
      fl   = (8.6*gamd + 1.35*(gam**(7./2.)))/(225. - 17.*gam + gamd )

        ! equation (4.5)
      if (abs(sqrt(gamd)) > 700.0) then
        qap=0.0
      else
        qap = fxy*(fl+ft)*exp(-gam)*(gamd**3)*(lamb**9)*3.0e21
      end if  
        q = qap*cvd
    end subroutine plasma

    subroutine pair(temp, rhom, qpai )
        ! pair.f
        ! modified from mkfit routine
        ! Itoh et al. ApJSS 102: 411 (1996)
        !
        !   input:      temp = temperature [K]
        !                       rhom = rho*Ye           [g cm**-3]
        ! output:       qpai = pair neutrino emissivity [ergs cm**-3 s**-1]
        !
        real(dp), intent(in) :: temp,rhom
        real(dp), intent(out) :: qpai
        real(dp) :: lam,xi,qp1,qp2,qp,fp1,fp2,fp,axi,g

        real(dp), parameter :: a0 =  6.002e19, a1 =  2.084e20, a2 =  1.872e21
        real(dp), parameter :: b1 =  9.383e-1, b2 = -4.141e-1, b3 =  5.829e-2, c =  5.5924e0
        real(dp), parameter :: bh1 =  1.2383, bh2 = -0.8141, bh3 =  0.0, ch =  4.9924
        real(dp), parameter :: coep = (cn_cv**2 + cn_ca**2) + cn_dn*(cn_cvp**2 + cn_cap**2)
        real(dp), parameter :: coem = (cn_cv**2 - cn_ca**2) + cn_dn*(cn_cvp**2 - cn_cap**2)

        ! equations (2.9)-(2.10)
        lam = temp/5.9302e9
        xi = ( rhom*1.0e-9 )**( 1.0/3.0 )/lam

        ! equation (2.6)
        qp1 = 10.74800*lam**2 + 0.3967*sqrt(lam) + 1.005
        qp2 = rhom/(7.692e7*lam**3 + 9.715e6*sqrt(lam))
        qp = (1.0/qp1)*(1.0 + qp2)**(-0.3)

        axi = a0+a1*xi+a2*xi**2

        ! equation (2.7)
        if (xi > 100.0) then
            fp1 = 0.0
        else if (temp < 1.0e10) then
            fp1 = axi*exp(-c*xi)
        else
            fp1 = axi*exp(-ch*xi)
        end if

        if (temp < 1.0e10) then
            fp2 = xi**3 + b1/lam + b2/( lam**2 ) + b3/( lam**3 )
        else
            fp2 = xi**3 + bh1/lam + bh2/( lam**2 ) + bh3/( lam**3 )
        end if
        fp = fp1/fp2

        ! equation (2.8)
        g = 1.0 - 13.04*lam**2 + 133.5*lam**4 +  1534.0*lam**6 + 918.6*lam**8

        ! equation (2.5)
        if ( lam < 7.0e-3 ) then
            qpai = 0.0
        else
            qpai = 0.5*coep*(1.0 + coem/coep*qp)*g*exp(-2.0/lam)*fp
        end if
    end subroutine pair

    subroutine photo(qpho,rhom,lamb)
        real(dp), intent(in) :: rhom,lamb
        real(dp), intent(out) :: qpho
        real(dp) :: qp,xi,fp
        ! equation (3.2)
        real(dp), parameter :: coe1 = 0.5*(( cn_cv**2 + cn_ca**2 ) + cn_dn*( cn_cvp**2 + cn_cap**2 ) )
        real(dp), parameter :: coe2 =  ((cn_cv**2 - cn_ca**2)  &
                & + cn_dn*(cn_cvp**2 - cn_cap**2))/((cn_cv**2 + cn_ca**2) + cn_dn*(cn_cvp**2 + cn_cap**2))

        qp = cqp(rhom,lamb)
        xi = ( rhom*1.e-9 )**( 1./3. )/lamb
        call cfp(fp,lamb,xi)
      qpho = coe1*(1.0 - coe2*qp)*rhom*lamb**5*fp
    
        contains
        function cqp(rhom,lamb) result(qp)
            real(dp), intent(in) :: rhom, lamb
            real(dp) :: qp
            qp = 0.666*(1.0 + 2.045*lamb)**(-2.066)     &
                        &   *(1.0 + rhom*(1.875e8*lamb + 1.653e8*lamb**2        &
                        &   + 8.499e8*lamb**3 - 1.604e8*lamb**4)**(-1))**(-1)
        end function cqp
        
    end subroutine photo
    
    subroutine cfp(fp,lamb,xi)
    ! cfp.f
    ! modified from cfp routine
    ! Itoh et al. ApJSS 102: 411 (1996)
    !
    ! input:    temp = temperature [K]
    !                   lamb = T/m_e
    !                   xi = normalized Fermi momentum
    !   returns: fp = fitting function (eq. [3.4]-[3.8])
    !
        real(dp),intent(in) :: lamb,xi
        real(dp),intent(out) :: fp
      real(dp) :: a(0:2), cc(0:2,0:6), d(0:2,1:5)
        real(dp) :: pjt,tau, c, temp
        real(dp), parameter :: pi = 3.141592654
        real(dp), parameter :: b(3) = (/ 6.290e-3, 7.483e-3, 3.061e-4 /) 
        integer :: i, j, k

        temp = lamb*5.9302e9
        if (temp < 1.0e7) then
            ! print *, 'WARNING: temperature < 1.0e7 K'
            fp = 0.0
            return
        else if (temp < 1.0e8) then
            include 'c7-8.par'
            include 'd7-8.par'
            c = 0.5654 + log10(temp/1.0e7)
            tau = log10(temp/1.0e7)
        else if (temp < 1.0e9) then
            include 'c8-9.par'
        include 'd8-9.par'
        c   =  1.5654
        tau = log10(temp/1.0e8)
      else
        include 'c9.par'
        include 'd9.par'
        c   =  1.5654
        tau = log10(temp/1.0e9)
      endif

        do j = 0, 2
            a(j)  = 0.5*cc(j,0) + 0.5*cc(j,6)*cos(10.0*pi*tau)
          do k = 1, 5
            pjt = pi*float(k)*tau
            a(j)= a(j) + cc(j,k)*cos(5.0/3.0*pjt)+d(j,k)*sin( 5.0/3.0*pjt )
            end do
        end do

        if ( c*xi > 100.0 ) then 
            fp = 0.d0
        endif

        fp = (a(0) + a(1)*xi + a(2)*xi**2)*exp(-c*xi) &
             &    / (xi**3 + b(1)/lamb+ b(2)/lamb**2 + b(3)/lamb**3) 
    end subroutine cfp
    
    subroutine pbf(kF,T,Tns,qpbf)
        ! pair-breaking and formation emissivity. Based on fmal. published by 
        ! Yakovlev & collaborators. 
        !
        ! Reduction factor follows recommendation of Steiner & Reddy (2009).
        use constants_def
        real(dp), intent(in) :: kF,T,Tns
        real(dp), intent(out) :: qpbf
        real(dp), parameter :: qpbf_pre = 1.17e21
        real(dp) :: mstar, xF, eF, R, tau, v, v2, F1, F2, F3, F, reduct
        
        tau = T/Tns
        qpbf = 0.0
        if (kF == 0.0 .or. tau > 1.0 .or. tau < 0.01 ) return
        xF = kF*hbarc_n/Mn_n
        eF = sqrt(1.0+xF**2)
        mstar = eF
        v = sqrt(1.0-tau)*(1.456-0.157/sqrt(tau)+1.764/tau)
        v2 = v**2
        F1 = v2*(0.602 + v2*(0.5942 + 0.288*v2))
        F2 = sqrt(0.5547+sqrt(0.4453**2+0.01130*v2))
        F3 = exp(-sqrt(4.0*v2+2.245**2)+2.245)
        F = F1*F2*F3
        reduct = (xF/eF)**2
        qpbf = qpbf_pre *(cn_dn+1)* mstar*xF*(T/1.0e9)**7 * F*reduct
    end subroutine pbf
    
    subroutine bremsstrahlung(rho,T,ionic,qbrems)
        use constants_def
        use nucchem_def, only: composition_info_type
        real(dp), intent(in) :: rho,T
        type(composition_info_type), intent(in) :: ionic
        real(dp), intent(out) :: qbrems
        real(dp), parameter :: threefourth = 3.0/4.0, GF = 1.436d-49, Cp2 = 1.675
        real(dp) :: qbrems_pre
        real(dp) :: re, kF, tau, eta, Zbareta, A, B, L
        
        qbrems_pre =  8.0*pi*GF**2*finestructure**2*Cp2*(boltzmann/hbar/clight)**6/567.0/hbar/amu
        re = (threefourth/pi*amu/rho/ionic%Ye)**onethird
        kF = (threepisquare*rho*ionic%Ye/amu)**onethird
        tau = 0.5*boltzmann*T/hbar/clight/kF
        eta = (ionic%A/ionic%Z)**onethird/re/cm_to_fm
        Zbareta = ionic%Z * eta
        A = 0.269 + 20.0*tau + 0.0168*ionic%Z + 0.00121*eta - 0.0356*Zbareta + 0.0137*ionic%Z2*tau + 1.54*Zbareta*tau
        B = 1.0 + 180.0*tau**2 + 0.483*tau*ionic%Z + 20.0*tau*Zbareta*eta + 4.31e-5*ionic%Z2
        L = A/B**threefourth
        qbrems = qbrems_pre*T**6*ionic%Z2*rho*L*(1.0-ionic%Yn)/ionic%A
    end subroutine bremsstrahlung
    
end module eval_crust_neutrino
