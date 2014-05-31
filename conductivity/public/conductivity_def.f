module conductivity_def
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

	type conductivity_components
		real :: total
		real :: ee
		real :: ei
		real :: eQ
		real :: sf
	end type conductivity_components

end module conductivity_def
