module conductivity_def
    use constants_def, only: dp

    ! used for mask array to control which components are included
    integer, parameter :: &
    &		icond_ee	= 1,	 &
    &		icond_ei	= 2,	 &
    &		icond_eQ	= 3,	 &
    &		icond_sf	= 4,     &
    &       icond_kap   = 5,	 &
    &		icond_nQ	= 6, 	 &
    &		icond_np	= 7
    integer, parameter :: num_conductivity_channels = 7
    logical, dimension(num_conductivity_channels), parameter :: cond_use_all =  &
    & [ .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE.]
    logical, dimension(num_conductivity_channels), parameter :: cond_use_only_conduction = &
    & [ .TRUE., .TRUE., .TRUE., .TRUE., .FALSE., .TRUE., .TRUE. ]
    logical, dimension(num_conductivity_channels), parameter :: cond_use_only_kap = &
    & [ .FALSE., .FALSE., .FALSE., .FALSE., .TRUE., .FALSE., .FALSE.]
    logical, dimension(num_conductivity_channels), parameter :: cond_exclude_sf = &
    & [ .TRUE., .TRUE., .TRUE., .FALSE., .TRUE., .TRUE., .TRUE. ]

    ! flag to control which ee scattering fmla. is used. 
    ! default: Shternin & Yakovlev (2006)
    integer, parameter :: icond_sy06 = 1
    ! Potekhin, Chabrier and Yakovlev (1997)
    integer, parameter :: icond_pcy = 2

    ! flag to control which eQ scattering fmla. is used.
    ! default: Potekhin (private comm.)
    integer, parameter :: icond_eQ_potekhin = 1
    ! Dany Page (private comm.)
    integer, parameter :: icond_eQ_page = 2

    type conductivity_components
        real(dp) :: total
        real(dp) :: ee
        real(dp) :: ei
        real(dp) :: eQ
        real(dp) :: sf
        real(dp) :: kap
        real(dp) :: nQ
        real(dp) :: np
        real(dp) :: electron_total
        real(dp) :: neutron_total
    end type conductivity_components

end module conductivity_def
