! extends MESA constants; checked against CODATA 2018

module constants_def
    use const_def
    
    implicit none

    ! mathematical constants
    real(dp), parameter :: twopi            = 2.0_dp*pi
    real(dp), parameter :: threepi          = 3.0_dp*pi
    real(dp), parameter :: fourpi           = pi4
    real(dp), parameter :: onethird         = one_third
    real(dp), parameter :: twothird         = two_thirds
    real(dp), parameter :: fivethird        = five_thirds
    real(dp), parameter :: seventhird       = 7.0_dp*onethird
    real(dp), parameter :: threepisquare    = (3.0_dp*pi**2)

    ! physical and astronomical constants, in cgs units
    real(dp), parameter :: boltzmann = kerg
    real(dp), parameter :: avogadro = avo
    real(dp), parameter :: clight2 = clight**2
    real(dp), parameter :: sigma_SB = boltz_sigma
    real(dp), parameter :: arad = crad
    real(dp), parameter :: si_charge = 1.602176634e-19_dp
    real(dp), parameter :: finestructure  = fine
    real(dp), parameter :: electroncharge2 = finestructure*hbar*clight
    real(dp), parameter :: electroncharge = sqrt(electroncharge2)
    real(dp), parameter :: Gnewton = standard_cgrav
    real(dp), parameter :: Mneutron = mn
    real(dp), parameter :: Mproton = mp
    real(dp), parameter :: Melectron = me
    real(dp), parameter :: Mmuon = 1.883531627e-25_dp
    real(dp), parameter :: theta_weak = 0.22290_dp  ! sin(theta_weak)**2
    real(dp), parameter :: Thomson = 8.0_dp*pi*onethird*(electroncharge2/Melectron/clight2)**2

    ! astronomical constants -- see http://ssd.jpl.nasa.gov/?constants
    real(dp), parameter :: julian_day = 86400_dp
    real(dp), parameter :: julian_year = 365.25_dp * julian_day
    real(dp), parameter :: GMsun = mu_sun

    ! some conversion factors
    real(dp), parameter :: fm_to_cm = 1.0e-13_dp
    real(dp), parameter :: ergs_to_mev = 1.0_dp/mev_to_ergs
    real(dp), parameter :: cm_to_fm = 1.0_dp/fm_to_cm
    
    ! atomic units (length in Bohr radii for hydrogen, energy in Rydberg)
    real(dp), parameter :: a_Bohr = hbar**2/electroncharge2/Melectron
    real(dp), parameter :: Rydberg = 0.5_dp*electroncharge2/a_Bohr
    
    ! nuclear units (c = 1, hbar = 197 MeV fm; energy is expressed in MeV, length in fm)
    real(dp), parameter :: hbarc_n = hbar*clight/mev_to_ergs/fm_to_cm
    real(dp), parameter :: mass_n = mev_to_ergs/clight2
    real(dp), parameter :: length_n = fm_to_cm
    real(dp), parameter :: density_n = 1.0_dp/fm_to_cm**3
    real(dp), parameter :: pressure_n= mev_to_ergs*density_n
    real(dp), parameter :: temperature_n = mev_to_ergs/boltzmann
    real(dp), parameter :: Mn_n = Mneutron*clight2*ergs_to_mev
    real(dp), parameter :: Mp_n = Mproton*clight2*ergs_to_mev
    real(dp), parameter :: Me_n = Melectron*clight2*ergs_to_mev
    real(dp), parameter :: amu_n = amu*clight2*ergs_to_mev
    
    ! gravitational units (G = c = 1; unit of mass is Msun)
    real(dp), parameter :: mass_g = Msun
    real(dp), parameter :: potential_g = clight2
    real(dp), parameter :: length_g = Gnewton*mass_g/potential_g
    real(dp), parameter :: density_g = mass_g/length_g**3
    real(dp), parameter :: pressure_g = density_g*potential_g
    real(dp), parameter :: time_g = length_g/clight

end module constants_def
