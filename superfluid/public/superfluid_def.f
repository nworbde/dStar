module superfluid_def
    use const_def, only: dp
    
    ! computes superfluid critical temperatures as a function of density
    ! these come courtesy of Dany Page

    ! types. these are pointers to locations in sf_tables
    integer, parameter :: undefined_gap = -1
    integer, parameter :: proton_1S0 = 1
    integer, parameter :: neutron_1S0 = 2
    integer, parameter :: neutron_3P2 = 3
    integer, parameter :: max_number_sf_types = 3
    
    type sf_table_type
        integer :: which_gap
		character(len=16) :: ref
        logical :: is_loaded
        integer :: nv
        real(dp) :: kF_min, kF_max
        real(dp), dimension(:), allocatable :: kF
        real(dp), dimension(:), pointer :: f=>null()        ! for storing the interpolation coefficients
    end type sf_table_type
    
    type (sf_table_type), dimension(max_number_sf_types), target :: sf_tables
    
    logical :: sf_is_initialized = .false.
    character(len=256), save :: sf_datadir
    
    ! you might want to scale the gaps by some factor; this is normally set to one
    real(dp), dimension(max_number_sf_types), save :: sf_scale
    
end module superfluid_def
