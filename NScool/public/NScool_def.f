module NScool_def
    use constants_def, only: dp
    use nucchem_def, only: composition_info_type
    
    ! storage for extra controls
    integer, parameter :: num_extra_real_controls = 32
    integer, parameter :: num_extra_integer_controls = 32
    integer, parameter :: num_extra_logical_controls = 32

    ! maximum number of sub-intervals
    integer, parameter :: max_number_epochs = 64

    ! interfaces for customizable routines
    abstract interface
    subroutine set_Qimp_interface(id,ierr)
        integer, intent(in) :: id
        integer, intent(out) :: ierr
    end subroutine set_Qimp_interface
    subroutine set_sf_interface(id,kp,kn,Tc)
        use constants_def
        use superfluid_def
        use superfluid_lib
        integer, intent(in) :: id   ! for accessing parameters
        real(dp), intent(in) :: kp, kn  ! proton, neutron wavevectors (fm**-1)
        real(dp), dimension(max_number_sf_types), intent(out) :: Tc ! (K)
    end subroutine set_sf_interface
    end interface

    type NScool_info
        integer :: id     ! id for this star
            
        include 'NScool_controls.inc'
        
        ! global information (NB. core here means interior to model)
        real(dp) :: Lsurf     ! emergent luminosity
        real(dp) :: dlnLsdlnT ! change surface luminosity with temperature out outer zone
        real(dp) :: Teff      ! surface effective temperature
        real(dp) :: Lcore     ! core luminosity
        real(dp) :: Tcore     ! core temperature
        real(dp) :: Psurf     ! surface pressure
        real(dp) :: Pcore     ! pressure at core boundary
        real(dp) :: Mcore     ! mass of core
        real(dp) :: Rcore     ! radius of core
        real(dp) :: ePhicore   ! potential of core
        real(dp) :: grav        ! gravity at surface
        
        real(dp) :: tsec      ! current value of time in seconds
        real(dp) :: dt        ! value of timestep just taken, in seconds
        integer :: model      ! counter that is incremented after each successful step
      
        ! information about the composition
        integer :: nisos  ! number of isotopes
        integer :: ncharged
!         integer, pointer, dimension(:) :: iso_ids ! id's for isotopes, computed by nucchem
        integer, pointer, dimension(:) :: charged_ids ! id's for charged isotopes, computed by nucchem
      
        ! administrative
        integer :: eos_handle
     
        character(len=256) :: base_profile_filename
        character(len=256) :: history_filename
        character(len=256) :: profile_manifest_filename
        
        ! zonal information
        integer :: nz     ! number of zones
        real(dp), pointer, dimension(:) :: dm      ! mass differences
!        real(dp), pointer, dimension(:,:) :: X     ! mass fracs (isotope, zone)
        real(dp), pointer, dimension(:,:) :: Yion  ! abundances of charged species
        real(dp), pointer, dimension(:) :: Xneut    ! neutron mass fraction
        type(composition_info_type), pointer, dimension(:) :: ionic ! composition info
        real(dp), pointer, dimension(:) :: P       ! pressure
        real(dp), pointer, dimension(:) :: lnT     ! ln(temperature)
        real(dp), pointer, dimension(:) :: T       ! temperature
        real(dp), pointer, dimension(:) :: ePhi     ! exp(potential)
        real(dp), pointer, dimension(:) :: e2Phi
        real(dp), pointer, dimension(:) :: eLambda ! redshift, 1+z
        real(dp), pointer, dimension(:) :: rho     ! density(P,T,X)
        real(dp), pointer, dimension(:) :: lnCp    ! ln(specific heat)
        real(dp), pointer, dimension(:) :: Cp      ! specific heat at const. pressure
        real(dp), pointer, dimension(:) :: dlnCp_dlnT ! derivative
        real(dp), pointer, dimension(:) :: lnGamma  ! ln(plasma Gamma)
        real(dp), pointer, dimension(:) :: Gamma    ! plasma Gamma: parameterizes strength of ion coulomb coupling
        real(dp), pointer, dimension(:) :: dlnGamma_dlnT    ! derivative
        real(dp), pointer, dimension(:) :: enu     ! specific neutrino emissivity
        real(dp), pointer, dimension(:) :: lnenu   ! ln(specific neutrino emissivity)
        real(dp), pointer, dimension(:) :: dlnenu_dlnT   ! derivative
        real(dp), pointer, dimension(:) :: enuc     ! heating rate (proportional to Mdot)

        ! facial information
        real(dp), pointer, dimension(:) :: m       ! mass
        real(dp), pointer, dimension(:) :: area    ! 4*pi*r**2
        real(dp), pointer, dimension(:) :: L       ! luminosity
        real(dp), pointer, dimension(:) :: dm_bar  ! interpolated mass difference
!         real(dp), pointer, dimension(:,:) :: X_bar  ! interpolated mass fracs
        real(dp), pointer, dimension(:,:) :: Yion_bar  ! interpolated abundances
        real(dp), pointer, dimension(:) :: Xneut_bar    ! interpolated neutron mass fraction
        type(composition_info_type), pointer, dimension(:) :: ionic_bar ! composition info
        real(dp), pointer, dimension(:) :: P_bar   ! interpolated pressure
        real(dp), pointer, dimension(:) :: lnT_bar ! ln(interpolated temperature)
        real(dp), pointer, dimension(:) :: T_bar   ! interpolated temperature
        real(dp), pointer, dimension(:) :: ePhi_bar ! exp(potential)
        real(dp), pointer, dimension(:) :: e2Phi_bar
        real(dp), pointer, dimension(:) :: eLambda_bar ! redshift, 1+z
        real(dp), pointer, dimension(:) :: rho_bar ! density(Pbar,Tbar,Xbar)
        real(dp), pointer, dimension(:) :: Kcond   ! thermal conductivity
        real(dp), pointer, dimension(:) :: lnK     ! ln(thermal conductivity)
        real(dp), pointer, dimension(:) :: dlnK_dlnT  ! derivative

        ! storage for interpolation of tabulated ln(enu), ln(Cp), ln(Kcond)
        integer :: n_tab                       ! number of table points
        real(dp), pointer, dimension(:) :: tab_lnT ! (n_tab) junction points for interpolation; same for all grid points
        real(dp), pointer, dimension(:,:) :: tab_lnenu  ! (4*n_tab, nz) coefficients for ln(enu)
        real(dp), pointer, dimension(:,:) :: tab_lnCp   ! (4*n_tab, nz) coefficients for ln(Cp)
        real(dp), pointer, dimension(:,:) :: tab_lnK    ! (4*n_tab, nz) coefficients for ln(Kcond)
        real(dp), pointer, dimension(:,:) :: tab_lnGamma  ! (4*n_tab, nz) coefficients for ln(plasma Gamma)

        ! storage for the lightcurve at selected points
        real(dp), pointer, dimension(:) :: t_monitor    ! (d, number_epochs)
        real(dp), pointer, dimension(:) :: Teff_monitor ! (K, number_epochs)

        logical :: in_use
        
        procedure(set_Qimp_interface), pointer, nopass ::  &
        & other_set_Qimp => null()
            
        procedure(set_sf_interface), pointer, nopass ::  &
        & other_sf_get_results => null()
        
    end type NScool_info

    integer, parameter :: max_NScool_handles = 10
    type(NScool_info), dimension(max_NScool_handles), target :: NScool_handles


contains
    subroutine get_NScool_info_ptr(id,s,ierr)
       integer, intent(in) :: id
       type(NScool_info), pointer, intent(out) :: s
       integer, intent(out) :: ierr
       if (id < 1 .or. id > max_NScool_handles) then
          ierr = -1
          return
       end if
       s => NScool_handles(id)
       ierr = 0
    end subroutine get_NScool_info_ptr
    
end module NScool_def
