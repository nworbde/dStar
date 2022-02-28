module dStar_crust_def
    use constants_def
    use nucchem_def, only: iso_name_length

    integer, parameter :: crust_default_number_table_points = 2048
    real(dp), parameter :: crust_default_lgPmin = 22.0
    real(dp), parameter :: crust_default_lgPmax = 33.5
    real(dp), parameter :: crust_default_temperature = 1.0e8_dp
    
    type crust_table_type
        logical :: is_loaded
        integer :: nv
        integer :: nisos
        character(len=iso_name_length), dimension(:), allocatable :: network
        real(dp) :: lgP_min, lgP_max
        real(dp) :: T   ! temperature at which EOS was computed
        real(dp), dimension(:), allocatable :: lgP
        ! storage of interpolation coefficients
        ! composition
        real(dp), dimension(:,:), pointer :: Y=>null()
        ! eos
        real(dp), dimension(:), pointer :: lgRho=>null() ! baryon density
        real(dp), dimension(:), pointer :: lgEps=>null() ! mass density
    end type crust_table_type

    type(crust_table_type), target :: crust_table
    logical, save :: crust_is_initialized = .FALSE.
    character(len=256), save :: crust_datadir
    
end module dStar_crust_def
