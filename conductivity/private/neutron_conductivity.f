module neutron_conductivity
    use conductivity_def
    implicit none
    
    ! for the integration over structure factors
    integer, parameter, private :: int_default_max_steps = 1000
    real(dp), parameter, private :: int_default_max_step_size = 0.0_dp
    real(dp), parameter, private :: int_default_starting_step = 1.0e-10_dp
    
contains

    function sPh(nn, nion, temperature, ionic)
        ! neutron superfluid phonon conductivity
        ! Implements formalism described in
        ! Aguilera et al. (2009), PRL 102: 091101
        !
        use math_lib
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

    function n_imp(nn, nion, temperature, ionic, ierr) result(nu)
        ! neutron-impurity scattering
        ! implements formalism described in Appendix of
        ! Deibel et al. (2017), ApJ 839: 95
        use exceptions_lib
        use nucchem_def, only: composition_info_type
        use constants_def
        use num_lib

        real(dp), intent(in) :: nn, nion, temperature
        type(composition_info_type), intent(in) :: ionic
        integer, intent(out) :: ierr
        real(dp) :: nu      ! neutron-impurity scattering frequency

        integer, parameter :: num_integration_variables = 1
        integer, parameter :: num_ipar = 0
        integer, parameter :: num_rpar = 6

        real(dp) :: ne      ! number density of electrons
        real(dp) :: kFn     ! neutron Fermi wavevector 
        real(dp) :: EFn     ! neutron Fermi energy
        real(dp) :: V       ! effective neutron-impurity potential 
        real(dp) :: Lambda_n_imp    ! Coulomb logarithm (Potekhin et al. 1999)
        real(dp) :: fac, sf_frac, Tc, xFn
        real(dp) :: gamma, vFn, v2, s, L1
        real(dp) :: pfn
        real(dp) :: R_a, mnstar
        real(dp) :: isospin, n_0, n_2, n_l, n_in    
        real(dp), dimension(:), pointer :: y
        integer ,dimension(:), pointer :: iwork => null()
        real(dp), dimension(:), pointer :: work => null()
        real(dp), dimension(num_integration_variables) :: rtol, atol
        integer, dimension(:), pointer :: ipar => null()
        real(dp), dimension(:), pointer :: rpar => null()
        real(dp) :: h, kn_end, kn_start, kn, kfi
        integer :: liwork, lwork, itol, lipar, lrpar, idid, lout, iout, npts
        character(len=128) :: msg
        type(alert) :: status=alert(scope='n_imp')

        ierr = 0
        allocate(y(num_integration_variables))

        kFn = (threepisquare*nn)**onethird
        kfi = (threepisquare*nn/ionic%Yn*(1.0-ionic%Yn)/ionic%A)**onethird

        call dop853_work_sizes( &
        &   num_integration_variables,num_integration_variables,liwork,lwork)
        allocate(iwork(liwork), work(lwork),ipar(num_ipar), rpar(num_rpar))
        iwork = 0
        iwork(5) = num_integration_variables
        work = 0.0
    
        itol = 0
        rtol = 1.0e-4
        atol = 1.0e-5
        iout = 0    ! solout is never called for output
        lout = error_unit
        lipar = num_ipar
        lrpar = num_rpar
                
        !calculate n_in and R_a using Sly4 model for now
        isospin = 1.0 - 2.0*ionic%Z/ionic%A
        n_0 = 0.1740/(fm_to_cm)**3
        n_2 = -0.0157/(fm_to_cm)**3
        n_l = n_0+n_2*isospin**2
        n_in = n_l/2.0*(1.0+0.9208*isospin)
        R_a = (3.0*(ionic%A-ionic%Z)/4.0/pi/n_in)**onethird    

        y(1) = 0.0
        kn_start = int_default_starting_step 
        kn_end = 2.0*kFn*R_a
        h = 0.0

        rpar(1) = temperature
        rpar(2) = ionic% Z
        rpar(3) = ionic% A
        rpar(4) = ionic% Yn
        rpar(5) = kfn
        rpar(6) = R_a

        call dop853(num_integration_variables, &
        &   structure_factor_impurity,kn_start,y,kn_end,h, &
        &   int_default_max_step_size,int_default_max_steps, &
        &   rtol,atol,itol, null_solout, iout, work, lwork, iwork, liwork,  &
        &   num_rpar, rpar, num_ipar, ipar, lout, idid)

        if (idid < 0) then
            write(msg,'(a,i0)') 'dop853 returned with error ',idid
            call status% report(msg)
            ierr = idid
        end if

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
        ! neutron-phonon scattering
        ! implements formalism described in Appendix of
        ! Deibel et al. (2017), ApJ 839: 95
        use exceptions_lib
        use nucchem_def, only: composition_info_type
        use constants_def
        use num_lib    

        real(dp), intent(in) :: nn, nion, temperature
        type(composition_info_type), intent(in) :: ionic
        integer, intent(out) :: ierr
        real(dp) :: nu              ! neutron-impurity scattering frequency

        integer, parameter :: num_integration_variables = 1
        integer, parameter :: num_ipar = 0
        integer, parameter :: num_rpar = 6
    
        real(dp) :: ne              ! electron number density
        real(dp) :: kFn             ! neutron Fermi wavevector 
        real(dp) :: EFn             ! neutron Fermi energy
        real(dp) :: V               ! effective neutron-impurity potential 
        real(dp) :: Lambda_n_phonon ! Coulomb logarithm (Potekhin et al. 1999)
        real(dp) :: fac, sf_frac, Tc, xFn
        real(dp) :: kfe, xr_e, eta, R_a, mnstar
        real(dp) :: qD, qi, ktf, beta, qs, wf, sf, L2, eta_n_0, G_k, Tp, eta_n
        real(dp) :: isospin, n_0, n_2, n_l, n_in
        real(dp), dimension(:), pointer :: y
        integer ,dimension(:), pointer :: iwork => null()
        real(dp), dimension(:), pointer :: work => null()
        real(dp), dimension(num_integration_variables) :: rtol, atol
        integer, dimension(:), pointer :: ipar => null()
        real(dp), dimension(:), pointer :: rpar => null()
        real(dp) :: h, kn_end, kn_start, kn, kfi
        integer :: liwork, lwork, itol, lipar, lrpar,idid, lout, iout, npts
        character(len=128) :: msg
        type(alert) :: status=alert(scope='n_phonon')

        ierr = 0
        allocate(y(num_integration_variables))

        kFn = (threepisquare*nn)**onethird
        kfi = (threepisquare*nn/ionic%Yn*(1.0-ionic%Yn)/ionic%A)**onethird

        call dop853_work_sizes( &
        &   num_integration_variables,num_integration_variables,liwork,lwork)
        allocate(iwork(liwork), work(lwork),ipar(num_ipar), rpar(num_rpar))
        iwork = 0
        iwork(5) = num_integration_variables
        work = 0.0
    
        itol = 0
        rtol = 1.0e-4
        atol = 1.0e-5
        iout = 0 ! solout is never called
        lout = error_unit
        lipar = num_ipar
        lrpar = num_rpar
              
        y(1) = 0.0

        !calculate n_in and R_a using Sly4 model for now
        isospin = 1.0 - 2.0*ionic%Z/ionic%A
        n_0 = 0.1740/(fm_to_cm)**3
        n_2 = -0.0157/(fm_to_cm)**3
        n_l = n_0+n_2*isospin**2
        n_in = n_l/2.0*(1.0+0.9208*isospin)
        R_a = (3.0*(ionic%A-ionic%Z)/4.0/pi/n_in)**onethird

        kn_start = int_default_starting_step
        kn_end = 2.0*kFn*R_a
        h = 0.0
        rpar(1) = temperature
        rpar(2) = ionic% Z
        rpar(3) = ionic% A
        rpar(4) = ionic% Yn
        rpar(5) = kfn
        rpar(6) = R_a

        call dop853(num_integration_variables, &
        &   structure_factor_phonon,kn_start,y,kn_end,h, &
        &   int_default_max_step_size,int_default_max_steps, &
        &   rtol,atol,itol, null_solout, iout, work, lwork, iwork, liwork,  &
        &   num_rpar, rpar, num_ipar, ipar, lout, idid)
    
        if (idid < 0) then
            write(msg,'(a,i0)') 'dop853 returned with error ',idid
            call status% report(msg)
            ierr = idid
        end if

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
        ! eq. (23), Deibel et al. (2017)
        use constants_def
        real(dp), intent(in) :: nn, n_in, Z, A
        real(dp) :: V
        V = hbar**2*(threepisquare*n_in)**(2.0/3.0)/2.0/Mneutron &
        &   *(1.0 - (nn/n_in)**(2.0/3.0))
    end function neutron_potential

    subroutine structure_factor_phonon(n,x,h,y,dy,lrpar,rpar,lipar,ipar,ierr)
        ! integrand for eq. (16), Deibel et al. (2017)
        use math_lib
        use constants_def

        integer, intent(in) :: n, lrpar, lipar
        real(dp), intent(in) ::  x,h
        real(dp), intent(inout) :: y(:)
        real(dp), intent(inout) :: dy(:)
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
        ! integrand for eq. (20), Deibel et al. (2017)
        use math_lib
        use constants_def

        integer, intent(in) :: n, lrpar, lipar
        real(dp), intent(in) ::  x,h
        real(dp), intent(inout) :: y(:)
        real(dp), intent(inout) :: dy(:)
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
    
    function core_neutron_conductivity(nn,np,mneff,mpeff,T,Tcs) result (Kn)
        use math_lib
        use constants_def
        use superfluid_def, only: max_number_sf_types, proton_1S0, neutron_1S0, neutron_3P2
        ! input nn, np: densities (fm**-3)
        ! T is temperature in K, Tcs are critical temperatures in K
        ! mneff = effective mass of neutron/neutron rest mass
        ! mpeff = effective mass of proton/proton rest mass
        real(dp), intent(in) :: nn,np,mneff,mpeff,T,Tcs(max_number_sf_types)
        real(dp) :: Kn
        real(dp) :: T8,tps,tns,tnt,istps,istns,vps,vns,vnt,yn,yp
        real(dp), parameter :: n0 = 0.16_dp  ! fm**-3
        real(dp) :: y, Sn1, Sn2, Sp1, Sp2
        real(dp) :: kFn, kFp, kFn2, kFp2
        real(dp) :: Kn1, Kn2, Kp1, Kp2, un1, un2, up1, up2
        real(dp) :: nu_nn, nu_np
        real(dp) :: mneffinv2, mpeffinv2

        mneffinv2 = 1.0/mneff**2
        mpeffinv2 = 1.0/mpeff**2

        tps = T/Tcs(proton_1S0)
        tns = T/Tcs(neutron_1S0)
        tnt = T/Tcs(neutron_3P2)

        istps = 1./sqrt(tps)
        istns = 1./sqrt(tns)
        vps = 0.0
        if(tps < 1.0) vps = sqrt(1.-tps)*(1.456+istps*(-0.157+istps*1.764))
        vns = 0.0
        if(tns < 1.0) vns = sqrt(1.-tns)*(1.456+istns*(-0.157+istns*1.764))
        vnt = 0.0
        if(tnt < 1.0) vnt = sqrt(1.-tnt)*(0.7893+1.188/tnt)

        yn = max(vns,vnt)
        yp = vps

        kFn = (threepisquare*nn)**onethird
        kFp = (threepisquare*np)**onethird
        kFn2 = kFn**2
        kFp2 = kFp**2

        Sn1 = 14.57/kFn**1.5 * (1.0-0.0788*kFn+0.0883*kFn2)/(1.0-0.1114*kFn)
        Sn2 = 7.880/kFn2 * (1.0-0.2241*kFn+0.2006*kFn2)/(1.0-0.1742*kFn)
        Sp1 = 0.8007*kFp/kFn**2 *  &
        &   (1.0+31.28*kFp-0.0004285*kFp**2 + 26.85*kFn + 0.08012*kFn**2)/ &
        &   (1.0-0.5898*kFn+0.2368*kFn**2+0.5838*kFp**2+0.884*kFn*kFp)
        Sp2 = 0.3830*kFp**4/kFn**5.5*(1.0+102.0*kFp+53.91*kFn)/ &
        &   (1.0-0.7087*kFn+0.2537*kFn**2+9.404*kFp**2-1.589*kFn*kFp)

        ! in-medium corrections
        !       un1 = kFn-1.665
        ! Kn1 = mneffinv2 * (0.4583+0.892*un1**2-0.5497*un1**3-0.06205*kFp+0.04022*kFp**2+0.2122*un1*kFp)
        ! un2 = kFn-1.556
        ! Kn2 = mneffinv2 * (0.4891+1.111*un2**2-0.2283*un2**3+0.01589*kFp-0.02099*kFp**2+0.2773*un2*kFp)
        !       up1 = kFn-2.126
        ! Kp1 = mpeffinv2 * (0.04377+1.1*up1**2+0.1180*up1**3+0.1626*kFp+0.3871*up1*kFp-0.2990*up1**4)
        ! up2 = kFn-2.116
        ! Kp2 = mpeffinv2 * (0.0001313+1.248*up2**2+0.2403*up2**3 &
        !    & +0.3257*kFp+0.5536*up2*kFp-0.3237*up2**4+0.09786*up2**2*kFp)

        Kn1 = 1.0_dp
        Kn2 = 1.0_dp
        Kp1 = 1.0_dp
        Kp2 = 1.0_dp

        ! put the whole thing together
        T8 = T*1.0e-8_dp
        nu_nn = 3.48e15*T8**2 *mneff**3 * &
        &      (Sn2*Kn2*Rn2(yn)+3.0*Sn1*Kn1*(Rn1(yn)-Rn2(yn)))
        nu_np = 3.48e15*T8**2*mneff*mpeff**2* &
        &      (Sp2*Kp2*Rp2(yn,yp)+0.5*Kp1*Sp1*(3.0*Rp1(yn,yp)-Rp2(yn,yp)))

        Kn = 7.2e23*T8*RC(yn)**2 /mneff * (1.0e15/(nu_nn+nu_np)) * nn/n0

    contains
        function RC(y)
            real(dp), intent(in) :: y
            real(dp) :: RC,y2
            y2 = y**2
            RC = (0.647+sqrt(0.353**2 + 0.109*y2))**1.5 &
            &       * exp(1.39-sqrt(1.39**2+y2))
        end function RC
        function Rn1(y)
            real(dp), intent(in) :: y
            real(dp) :: Rn1, y2
            y2 = y**2
            Rn1 = twothird*(0.9468+sqrt(0.0532**2+0.5346*y2))**3 &
            & * exp(0.377-sqrt(0.377**2 + 4.0*y2)) &
            & + onethird*(1.0+1.351*y2)**2 &
            & * exp(0.169-sqrt(0.169**2+9.0*y2))
        end function Rn1
        function Rn2(y)
            real(dp), intent(in) :: y
            real(dp) :: Rn2,y2
            y2 = y**2
            Rn2 = 0.5*(0.6242+sqrt(0.3758**2 + 0.07198*y2))**2  &
            & * exp(3.6724 - sqrt(3.6724**2+4.0*y2))  &
            & + 0.5*(1.0+0.01211*y2)**9 &
            & * exp(7.5351-sqrt(7.5351**2+9.0*y2))
        end function Rn2
        function Rp1(yn,yp)
            real(dp), intent(in) :: yn,yp
            real(dp) :: Rp1
            real(dp) :: yn2,yp2,uplus, uminus, un, up, yminus, yplus
            yn2 = yn**2
            yp2 = yp**2
            if (yn2 > 0.0 .and. yp2 == 0.0) then
                Rp1 = (0.4459 + sqrt(0.5541**2 + 0.03016*yn2))**2 &
                & * exp(2.1178-sqrt(2.1178**2 + yn2))
            else if (yn2 == 0.0 .and. yp2 > 0.0) then
                Rp1 = 0.5*(0.3695+sqrt(0.6305**2+0.01064*yp2)) &
                & * exp(2.4451-sqrt(2.4451**2+yp2)) &
                & + 0.5*(1.0+0.1917*yp2)**1.4 &
                & * exp(4.6627-sqrt(4.6627**2+4.0*yp2))
            else if (yn2 > 0.0 .and. yp2 > 0.0) then
                yminus = min(yn, yp)
                yplus = max(yn, yp)
                uplus = sqrt(yplus**2+1.485**2)-1.485
                uminus = sqrt(yminus**2+1.485**2)-1.485
                un = sqrt(yn2+1.485**2)-1.485
                up = sqrt(yp2+1.485**2)-1.485
                Rp1 = exp(-uplus-uminus)*(0.7751+0.4823*un+0.1124*up &
                & + 0.04991*un**2+0.08513*un*up+0.01284*un**2*up) &
                & + exp(-2.0*uplus) *(0.2249+0.3539*uplus &
                &       -0.2189*uminus-0.6069*un*uminus + 0.7362*up*uplus)
            else
                Rp1 = 1.0
            end if
        end function Rp1
        function Rp2(yn,yp)
            real(dp), intent(in) :: yn,yp
            real(dp) :: Rp2, yn2, yp2
            real(dp) :: uplus,uminus,un,up,yminus,yplus
            yn2 = yn**2
            yp2 = yp**2
            if (yn2>0.0 .and. yp2 == 0.0) then
                Rp2 = (0.801+sqrt(0.199**2+0.04645*yn2))**2  &
                & * exp(2.3569-sqrt(2.3569**2+yn2))
            else if (yn2 == 0.0 .and. yp2 > 0.0) then
                Rp2 = 0.0436*(sqrt(4.345**2+19.55*yp2)-3.345) &
                & * exp(2.0247-sqrt(2.0247**2 + yp2)) &
                & + 0.0654*exp(8.992-sqrt(8.992**2+1.5*yp2)) &
                & + 0.891*exp(9.627-sqrt(9.627**2+9.0*yp2))
            else if (yn2 > 0.0 .and. yp2 > 0.0) then
                yplus = max(yn,yp)
                yminus = min(yn,yp)
                uplus = sqrt(yplus**2+1.761**2)-1.761
                uminus = sqrt(yminus**2+1.761**2)-1.761
                un  = sqrt(yn2+1.761**2)-1.761
                up = sqrt(yp2+1.761**2)-1.761
                Rp2 = exp(-uplus-uminus)*(1.1032+0.8645*un+0.2042*up &
                & + 0.07937*un**2+0.1451*un*up+0.01333*un**2*up) &
                & + exp(-2.0*uplus)*(-0.1032-0.2340*uplus+0.06152*un*uplus  &
                & + 0.7533*un*uminus-1.007*up*uplus)
            else
                Rp2 = 1.0
            end if
        end function Rp2
    end function core_neutron_conductivity
    
end module neutron_conductivity
