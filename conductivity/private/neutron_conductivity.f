module neutron_conductivity
    use conductivity_def
    implicit none
    
    ! for the integration over structure factors
    integer, parameter :: int_default_max_steps = 1000
    real(dp), parameter :: int_default_max_step_size = 0.0_dp
    real(dp), parameter :: int_default_starting_step = 1.0e-10_dp
    
contains

    function sPh(nn, nion, temperature, ionic)
        ! neutron superfluid phonon conductivity
        ! Implements formalism described in
        ! Aguilera et al. (2009), PRL 102: 091101
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

    function n_imp(nn, nion, temperature, ionic, ierr) result(nu)
        ! neutron-impurity scattering
        ! implements formalism described in Appendix of
        ! Deibel et al. (2017), ApJ 839: 95
        use, intrinsic :: iso_fortran_env, only: error_unit
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
            write(error_unit,*) 'n_imp: dop853 returned with error',idid
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
        use, intrinsic :: iso_fortran_env, only: error_unit
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
            write(error_unit,*) 'n_phonon: dop853 returned with error',idid
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
end module neutron_conductivity
