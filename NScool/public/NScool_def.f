module NScool_def

    integer, parameter :: num_extra_real_controls = 32
    integer, parameter :: num_extra_integer_controls = 32
    integer, parameter :: num_extra_logical_controls = 32

    type NScool_info
        integer :: id     ! id for this star
            
        include 'controls.inc'
        
        ! global information (NB. core here means interior to model)
        real(dp) :: Lsurf     ! emergent luminosity
        real(dp) :: dlnLsdlnT ! change surface luminosity with temperature out outer zone
        real(dp) :: Teff      ! surface effective temperature
        real(dp) :: Lcore     ! core luminosity
        real(dp) :: Psurf     ! surface pressure
        real(dp) :: Pcore     ! pressure at core boundary
        real(dp) :: tsec      ! current value of time in seconds
        real(dp) :: dt        ! value of timestep just taken, in seconds
        integer :: model      ! counter that is incremented after each successful step
      
        ! information about the composition
        integer :: nisos  ! number of isotopes
        integer :: ncharged
        integer, allocatable, dimension(:) :: iso_ids ! id's for isotopes, computed by nucchem
        integer, allocatable, dimension(:) :: charged_ids ! id's for charged isotopes, computed by nucchem
      
        ! administrative
        integer :: eos_handle
     
        character(len=256) :: base_profile_filename
        character(len=256) :: history_filename
      
        ! zonal information
        integer :: nz     ! number of zones
        real(dp), allocatable, dimension(:) :: dm      ! mass differences
        real(dp), allocatable, dimension(:,:) :: X     ! mass fracs (isotope, zone)
        real(dp), allocatable, dimension(:,:) :: Yion  ! abundances of charged species
        type(composition_info_type), allocatable, dimension(:) :: ionic ! composition info
        real(dp), allocatable, dimension(:) :: P       ! pressure
        real(dp), allocatable, dimension(:) :: lnT     ! ln(temperature)
        real(dp), allocatable, dimension(:) :: T       ! temperature
        real(dp), allocatable, dimension(:) :: phi     ! potential
        real(dp), allocatable, dimension(:) :: eLambda ! redshift, 1+z
        real(dp), allocatable, dimension(:) :: rho     ! density(P,T,X)
        real(dp), allocatable, dimension(:) :: lnCp    ! ln(specific heat)
        real(dp), allocatable, dimension(:) :: Cp      ! specific heat at const. pressure
        real(dp), allocatable, dimension(:) :: dlnCp_dlnT ! derivative
        real(dp), allocatable, dimension(:) :: enu     ! specific neutrino emissivity
        real(dp), allocatable, dimension(:) :: lnenu   ! ln(specific neutrino emissivity)
        real(dp), allocatable, dimension(:) :: dlnenu_dlnT   ! derivative

        ! facial information
        real(dp), allocatable, dimension(:) :: m       ! mass
        real(dp), allocatable, dimension(:) :: L       ! luminosity
        real(dp), allocatable, dimension(:) :: dm_bar  ! interpolated mass difference
        real(dp), allocatable, dimension(:,:) :: X_bar  ! interpolated mass fracs
        real(dp), allocatable, dimension(:,:) :: Yion_bar  ! interpolated abundances
        type(composition_info_type), allocatable, dimension(:) :: ionic_bar ! composition info
        real(dp), allocatable, dimension(:) :: P_bar   ! interpolated pressure
        real(dp), allocatable, dimension(:) :: lnT_bar ! ln(interpolated temperature)
        real(dp), allocatable, dimension(:) :: T_bar   ! interpolated temperature
        real(dp), allocatable, dimension(:) :: phi_bar ! potential
        real(dp), allocatable, dimension(:) :: eLambda_bar ! redshift, 1+z
        real(dp), allocatable, dimension(:) :: rho_bar ! density(Pbar,Tbar,Xbar)
        real(dp), allocatable, dimension(:) :: Kcond   ! thermal conductivity
        real(dp), allocatable, dimension(:) :: lnK     ! ln(thermal conductivity)
        real(dp), allocatable, dimension(:) :: dlnK_dlnT  ! derivative

        ! storage for interpolation of tabulated ln(enu), ln(Cp), ln(Kcond)
        integer :: n_tab                       ! number of table points
        real(dp), pointer, dimension(:) :: tab_lnT ! (n_tab) junction points for interpolation; same for all grid points
        real(dp), pointer, dimension(:,:) :: lnenu_tab  ! (4*n_tab, nz) coefficients for ln(enu)
        real(dp), pointer, dimension(:,:) :: lnCp_tab  ! (4*n_tab, nz) coefficients for ln(Cp)
        real(dp), pointer, dimension(:,:) :: lnK_tab   ! (4*n_tab, nz) coefficients for ln(Kcond)

        logical :: in_use
        
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
