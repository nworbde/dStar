module dStar_crust_def
    use constants_def

    integer, parameter :: crust_default_number_table_points = 2048
    real(dp), parameter :: crust_default_lgPmin = 26.5
    real(dp), parameter :: crust_default_lgPmax = 33.5
    
    type crust_table_type
        logical :: is_loaded
        integer :: nv
        real(dp) :: lgP_min, lgP_max
        real(dp), dimension(:), allocatable :: lgP
        ! storage of interpolation coefficients
        real(dp), dimension(:), pointer :: lgRho=>null()
        real(dp), dimension(:), pointer :: lgEps=>null()
    end type crust_table_type

    type(crust_table_type), target :: crust_table
    logical :: crust_is_initialized = .FALSE.
    character(len=256) :: crust_datadir
    
end module dStar_crust_def
