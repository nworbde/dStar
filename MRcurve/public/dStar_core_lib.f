module dStar_core_lib
	use dStar_core_def

contains
	
	subroutine dstar_core_startup(datadir,ierr)
		use, intrinsic :: iso_fortran_env, only: error_unit
		character(len=*), intent(in) :: datadir
		integer, intent(out) :: ierr
		if (dStar_core_is_initialized) then
			write (error_unit, *) &
			& 'dstar_core_startup called on an already initialized module'
			return
		end if
		core_datadir = trim(datadir)//'/skyrme_data'
		dStar_core_is_initialized = .TRUE.
		ierr = 0
	end subroutine dStar_core_startup
	
	subroutine dstar_core_shutdown()
		dStar_core_is_initialized = .FALSE.
		core_datadir = ''
	end subroutine dstar_core_shutdown
	
end module dStar_core_lib
