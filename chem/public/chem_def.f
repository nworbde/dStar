module chem_def
    use constants_def, only: dp
    use utils_def, only: integer_dict
    use netJina_def, only: nuclib_data
    implicit none

    integer, parameter :: nuclide_not_found = -1
    
	! for storage of composition information
	type composition_info_type
		real(dp) :: A
		real(dp) :: Z
		real(dp) :: Z53
		real(dp) :: Z2
		real(dp) :: Z73
		real(dp) :: Z52
		real(dp) :: ZZ1_32
		real(dp) :: Z2XoA2
		real(dp) :: Ye
		real(dp) :: Yn
		real(dp) :: Q
	end type composition_info_type

    type(integer_dict), pointer :: nuclide_dict=>null()
    type(nuclib_data) :: nuclib
    logical, save :: chem_is_initialized
    
end module chem_def
