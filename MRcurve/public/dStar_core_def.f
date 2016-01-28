module dStar_core_def
    use constants_def
	real(dp), parameter :: rho_saturation = 0.16	! fm
    integer, parameter :: core_default_number_table_points = 2048
    real(dp), parameter :: core_default_Rhomin = 0.25*rho_saturation
    real(dp), parameter :: core_default_Rhomax = 12.0*rho_saturation
    
    type core_table_type
        logical :: is_loaded
        integer :: nv
        real(dp) :: lgP_min, lgP_max
        real(dp), dimension(:), allocatable :: lgP
        ! storage of interpolation coefficients
        real(dp), dimension(:), pointer :: lgRho=>null()
        real(dp), dimension(:), pointer :: lgEps=>null()
		real(dp), dimension(:), pointer :: Xprot=>null()
    end type core_table_type

    type(core_table_type), target :: core_table
    
	integer, parameter :: skyrme_matter = 1, skyrme_neutron = 2, skyrme_sym = 3

	type skyrme
		integer :: model
		real(dp) :: gamma
		real(dp) :: a, b, c, d, e
	end type skyrme

	type(skyrme), target :: skyrme_eos(3)
	
	character(len=256), save :: core_datadir
	logical :: dStar_core_is_initialized = .FALSE.
	
end module dStar_core_def
