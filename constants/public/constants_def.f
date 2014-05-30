! extends MESA constants and updates them to CODATA 2010 values

module constants_def
    use const_def
    
    implicit none

    ! mathematical constants
    real(dp), parameter :: twopi            = 2.0_dp*pi
    real(dp), parameter :: threepi          = 3.0_dp*pi
    real(dp), parameter :: fourpi           = 4.0_dp*pi
    real(dp), parameter :: onethird         = one_third
    real(dp), parameter :: twothird         = two_thirds
    real(dp), parameter :: fivethird        = 5.0_dp*onethird
    real(dp), parameter :: seventhird       = 7.0_dp*onethird
    real(dp), parameter :: threepisquare    = (3.0_dp*pi**2)

    ! physical and astronomical constants, in cgs units
    real(dp), save :: boltzmann        ! = 1.3806488e-16_dp
    real(dp), save :: avogadro         ! = 6.02214129e23_dp
    real(dp), save :: clight2                  ! = clight**2
    real(dp), save :: sigma_SB                 ! = 5.670400e-5
    real(dp), save :: arad                         ! = sigma_SB*4.0/clight
    real(dp), save :: finestructure        ! = 7.2973525698e-3_dp
    real(dp), save :: Gnewton                  ! = 6.67384e-8_dp
    real(dp), save :: Mneutron                 ! = 1.674927351e-24_dp
    real(dp), save :: Mproton                  ! = 1.672621777e-24_dp
    real(dp), save :: si_charge                ! = 1.602176565e-19_dp
    real(dp), save :: electroncharge       ! = si_charge*clight*0.1_dp
    real(dp), save :: Melectron                ! = 9.10938291e-28_dp
    real(dp), save :: Mmuon                        ! = 1.883531475e-25_dp
    real(dp), save :: theta_weak               ! = 0.2223_dp

    ! astronomical constants -- see http://ssd.jpl.nasa.gov/?constants
    real(dp), save :: julian_day               ! = 86400_dp
    real(dp), save :: julian_year          ! = 365.25_dp * julian_day
    real(dp), save :: GMsun                        ! = 1.32712440018e26_dp

    ! some conversion factors
    real(dp), save :: fm_to_cm         ! = 1.0e-13_dp
    real(dp), save :: ergs_to_mev          ! = 1.0_dp/mev_to_ergs
    real(dp), save :: cm_to_fm                 ! = 1.0_dp/fm_to_cm
    
    ! atomic units (length in Bohr radii for hydrogen, energy in Rydberg)
    real(dp), save :: a_Bohr                       ! = hbar**2/Melectron/electroncharge**2
    real(dp), save :: Rydberg                  ! = 0.5_dp*electroncharge**2/a_Bohr
    
    ! nuclear units (c = 1, hbar = 197 MeV fm; energy is expressed in MeV, length in fm)
    real(dp), save :: hbarc_n              ! = hbar*clight/mev_to_ergs/fm_to_cm
    real(dp), save :: mass_n                       ! = mev_to_ergs/clight2
    real(dp), save :: length_n                 ! = fm_to_cm
    real(dp), save :: density_n                ! = 1.0_dp/fm_to_cm**3
    real(dp), save :: pressure_n               ! = mev_to_ergs*density_n
    real(dp), save :: temperature_n        ! = mev_to_ergs/boltzmann
    real(dp), save :: Mn_n                         ! = Mneutron*clight2*ergs_to_mev
    real(dp), save :: Mp_n                         ! = Mproton*clight2*ergs_to_mev
    real(dp), save :: Me_n                         ! = Melectron*clight2*ergs_to_mev
    real(dp), save :: amu_n                        ! = amu*clight2*ergs_to_mev
    
    ! gravitational units (G = c = 1; unit of mass is Msun)
    real(dp), save :: mass_g                       ! = Msun
    real(dp), save :: potential_g          ! = clight2
    real(dp), save :: length_g                 ! = Gnewton*mass_g/potential_g
    real(dp), save :: density_g                ! = mass_g/length_g**3
    real(dp), save :: pressure_g               ! = density_g*potential_g
    real(dp), save :: time_g                       ! = length_g/clight


contains

    subroutine initialize_constants

        boltzmann           = 1.3806488e-16_dp
        avogadro            = 6.02214129e23_dp
        clight              = 2.99792458e10_dp
        clight2             = clight**2
        sigma_SB            = 5.670400e-5
        arad                = sigma_SB*4.0/clight
        amu                 = 1.660538921e-24_dp
        hbar                = 1.054571726e-27_dp
        finestructure       = 7.2973525698e-3_dp
        Gnewton             = 6.67384e-8_dp
        Mneutron            = 1.674927351e-24_dp
        Mproton             = 1.672621777e-24_dp
        si_charge           = 1.602176565e-19_dp
        electroncharge      = si_charge*clight*0.1_dp
        Melectron           = 9.10938291e-28_dp
        Mmuon               = 1.883531475e-25_dp
        theta_weak          = 0.2223_dp

        julian_day          = 86400_dp
        julian_year         = 365.25_dp * julian_day
        GMsun               = 1.32712440018e26_dp
        Msun                = GMsun/Gnewton

        mev_to_ergs         = 1.602176565e-6_dp
        fm_to_cm            = 1.0e-13_dp
        ergs_to_mev         = 1.0_dp/mev_to_ergs
        cm_to_fm            = 1.0_dp/fm_to_cm

        a_Bohr              = hbar**2/Melectron/electroncharge**2
        Rydberg             = 0.5_dp*electroncharge**2/a_Bohr

        hbarc_n             = hbar*clight/mev_to_ergs/fm_to_cm
        mass_n              = mev_to_ergs/clight2
        length_n            = fm_to_cm
        density_n           = 1.0_dp/fm_to_cm**3
        pressure_n          = mev_to_ergs*density_n
        temperature_n       = mev_to_ergs/boltzmann
        Mn_n                = Mneutron*clight2*ergs_to_mev
        Mp_n                = Mproton*clight2*ergs_to_mev
        Me_n                = Melectron*clight2*ergs_to_mev
        amu_n               = amu*clight2*ergs_to_mev

        mass_g              = Msun
        potential_g         = clight2
        length_g            = Gnewton*mass_g/potential_g
        density_g           = mass_g/length_g**3
        pressure_g          = density_g*potential_g
        time_g              = length_g/clight

        standard_cgrav = Gnewton
        ! gravitational constant (g^-1 cm^3 s^-2)
        planck_h = hbar*2.0*pi 
        ! Planck's constant (erg s)
        qe = electroncharge
        ! electron charge (esu == (g cm^3 s^-2)^(1/2))
        avo = avogadro 
        ! Avogadro's constant (mole^-1)
        kerg = boltzmann 
        ! Boltzmann's constant (erg K^-1)
        boltzm = kerg
        cgas = boltzm*avo !ideal gas constant; erg/K
        kev = 8.617385d-5 
        ! converts temp to ev (ev K^-1)

        mn = Mneutron! neutron mass (g)
        mp = Mproton ! proton mass (g)
        me = Melectron ! (was 9.1093897d-28) electron mass (g)
        rbohr = a_Bohr ! Bohr radius (cm)
        fine = finestructure ! fine structure constant
        ev2erg = 1.602176487d-12 ! electron volt (erg)
        hion = Rydberg/ev2erg ! hydrogen ionization energy (eV)
        mev_amu = mev_to_ergs/amu
        Qconv = mev_to_ergs*avo

        boltz_sigma = sigma_SB
        ! boltzmann's sigma = a*c/4 (erg cm^-2 K^-4 s^-1)
        crad = arad
        ! = radiation density constant, a (erg cm^-3 K^-4); Prad = crad * T^4 / 3
        ! approx = 7.5657e-15

        ssol = boltz_sigma
        asol = crad
        weinlam = planck_h*clight/(kerg * 4.965114232d0)
        weinfre = 2.821439372d0*kerg/planck_h
        rhonuc = 2.342d14 ! density of nucleus (g cm^3)

        ! solar age, L, and R values from Bahcall et al, ApJ 618 (2005) 1049-1056.
        msol = Msun  ! solar mass (g)  <<< gravitational mass, not baryonic
        rsol = 6.9598d10 ! solar radius (cm)
        lsol = 3.8418d33  ! solar luminosity (erg s^-1)
        agesol = 4.57d9  ! solar age (years)
        Rsun = rsol
        Lsun = lsol
        Msun33 = msol*1d-33
        Rsun11 = rsol*1d-11
        Lsun33 = lsol*1d-33
        ly = 9.460528d17 ! light year (cm)
        pc = 3.261633d0 * ly ! parsec (cm)
        secyer = 3.1558149984d7 ! seconds per year

        m_earth = 5.9764d27 ! earth mass (g)
        ! = 3.004424e-6 Msun
        r_earth = 6.37d8 ! earth radius (cm)
        au = 1.495978921d13 ! astronomical unit (cm)

        m_jupiter = 1.8986d30 ! jupiter mass (g)
        ! = 0.954454d-3 Msun
        r_jupiter = 6.9911d9 ! jupiter mean radius (cm)
        semimajor_axis_jupiter = 7.7857d13 ! jupiter semimajor axis (cm)
    
    end subroutine initialize_constants

end module constants_def
