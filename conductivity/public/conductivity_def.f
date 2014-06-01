module conductivity_def
    use constants_def, only: dp
    
	! used for mask array to control which components are included
	integer, parameter :: &
			&		icond_ee	= 1,	 &
			&		icond_ei	= 2,	 &
			&		icond_eQ	= 3,	 &
			&		icond_sf	= 4
	integer, parameter :: num_conductivity_channels = 4
	logical, dimension(num_conductivity_channels), parameter :: cond_use_all =  &
				& [ .TRUE., .TRUE., .TRUE., .TRUE.]
	
	! flag to control which ee scattering fmla. is used. default is Shternin & Yakovlev '06
	integer, parameter :: icond_sy06 = 1, icond_pcy = 2

	! flag to control which eQ scattering fmla. is used. default is Potekhin
	integer, parameter :: icond_eQ_potekhin = 1, icond_eQ_page = 2

	type conductivity_components
		real(dp) :: total
		real(dp) :: ee
		real(dp) :: ei
		real(dp) :: eQ
		real(dp) :: sf
	end type conductivity_components

end module conductivity_def
