module neutrino_def
    use constants_def, only : dp
    ! routines for evaluating neutrino emissivity from uniform nuclear matter.  
    ! the rates use analytical fmla. published by Yakovlev
    ! and collaborators.

    ! the routines take a logical array to control which channels are included 
    integer, parameter ::  &
        & icore_nu_brem     = 1,  &
        & icore_nu_mUrca    = 2,    &
        & icore_nu_dUrca    = 3,    &
        & icore_nu_PBF      = 4
    integer, parameter :: num_core_nu_channels = 4

    logical, dimension(num_core_nu_channels), parameter :: &
    &    core_nu_minimal_cooling = [ .TRUE.,.TRUE.,.FALSE.,.TRUE. ]

    logical, dimension(num_core_nu_channels), parameter :: &
    &    core_nu_enhanced_cooling = [ .TRUE.,.TRUE.,.TRUE.,.TRUE. ]
    
    type core_neutrino_emissivity_channels
        real(dp) :: total   ! total neutrino emissivity (units?)
        real(dp) :: brem    ! nucleon-nucleon bremsstrahlung
        real(dp) :: brem_np ! neutron-proton bremsstrahlung
        real(dp) :: brem_pp ! proton-proton bremsstrahlung
        real(dp) :: brem_nn ! neutron-neutron bremsstrahlung
        real(dp) :: mUrca   ! mod. Urca
        real(dp) :: mUrca_p ! mod. Urca, proton branch
        real(dp) :: mUrca_n ! mod. Urca, neutron branch
        real(dp) :: dUrca   ! dir. Urca
        real(dp) :: PBF     ! PBF total
        real(dp) :: PBF_n   ! PBF, neutron
        real(dp) :: PBF_p   ! PBF, proton
    end type core_neutrino_emissivity_channels

    integer, parameter ::   &
        &   icrust_nu_pair      = 1,  &
        &   icrust_nu_photo     = 2,    &
        &   icrust_nu_plasma    = 3,    &
        &   icrust_nu_brems     = 4,    &
        &   icrust_nu_pbf       = 5
    integer, parameter :: num_crust_nu_channels = 5

    logical, dimension(num_crust_nu_channels), parameter :: &
    &    crust_nu_full_cooling = [ .TRUE., .TRUE., .TRUE., .TRUE., .TRUE. ] 

    type crust_neutrino_emissivity_channels
        real(dp) :: total
        real(dp) :: pair
        real(dp) :: photo
        real(dp) :: plasma
        real(dp) :: brems
        real(dp) :: pbf
    end type crust_neutrino_emissivity_channels

end module neutrino_def
