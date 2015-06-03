module dStar_atm_def
    use constants_def
    
	integer, parameter :: atm_default_number_table_points = 256
	real(dp), parameter :: atm_default_lgTbmin = 7.0
	real(dp), parameter :: atm_default_lgTbmax = 9.5
    real(dp), parameter :: default_lgTeff_min = 5.5
    real(dp), parameter :: default_lgTeff_max = 7.0
    

	type atm_table_type
		logical :: is_loaded
		integer :: nv
		real(dp) :: lgTb_min, lgTb_max
		real(dp), dimension(:), allocatable :: lgTb
        ! storage of interpolation coefficients
		real(dp), dimension(:), pointer :: lgTeff=>null()
		real(dp), dimension(:), pointer :: lgflux=>null()
	end type atm_table_type

	type(atm_table_type), target :: atm_table
	logical :: atm_is_initialized = .FALSE.
	character(len=256) :: atm_datadir

end module dStar_atm_def
