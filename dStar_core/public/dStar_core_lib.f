module dStar_core_lib
	use dStar_core_def

contains
	
	subroutine dstar_core_startup(datadir,ierr)
		use, intrinsic :: iso_fortran_env, only: error_unit
		character(len=*), intent(in) :: datadir
		integer, intent(out) :: ierr
		if (dStar_core_is_initialized) then
			ierr = 1
			write (error_unit, *) 'dstar_core_startup: ', &
			&	'package already initialized'
			return
		end if
		core_datadir = trim(datadir)//'/core_data'
		dStar_core_is_initialized = .TRUE.
		ierr = 0
	end subroutine dStar_core_startup
	
	subroutine dstar_core_shutdown()
		use dStar_core_mod, only : do_free_core_table
		type(core_table_type), pointer :: tab
		tab => core_table
		call do_free_core_table(tab)
		dStar_core_is_initialized = .FALSE.
	end subroutine dstar_core_shutdown
	
	subroutine dStar_core_free_table()
		use dStar_core_mod, only: do_free_core_table
		type(core_table_type), pointer :: tab
		tab => core_table
		call do_free_core_table(tab)
	end subroutine dStar_core_free_table
	
	subroutine dStar_core_load_table(eos_handle,ierr)
		use iso_fortran_env, only : error_unit
		use dStar_core_mod, only : do_load_core_table
		integer, intent(in) :: eos_handle
		integer, intent(out) :: ierr
		
		call do_load_core_table(eos_handle,ierr)
		if (ierr /= 0) then
			write(error_unit,'(a,i3)') 'dStar_core_load_table: ierr = ',ierr
		end if
	end subroutine dStar_core_load_table
	
	subroutine dStar_core_get_results(lgP,lgRho,dlgRho,lgEps,dlgEps,Xprot,dXprot,ierr)
		use iso_fortran_env, only : error_unit
		use interp_1d_lib, only: interp_value_and_slope
		real(dp), intent(in) :: lgP
		real(dp), intent(out) :: lgRho,dlgRho,lgEps,dlgEps,Xprot,dXprot
		integer, intent(out) :: ierr
		real(dp) :: lgPr
		type(core_table_type), pointer :: tab
		character(len=*), parameter :: routine_name='dStar_core_get_results'
		
		tab => core_table
		if (.not.tab% is_loaded) then
			ierr = -9
			write(error_unit,*) routine_name//': table is not loaded'
			return
		end if
		
		! clip lgP to table
		lgPr = max(lgP,tab% lgP_min)
		lgPr = min(lgP,tab% lgP_max)
		call interp_value_and_slope(tab% lgP, tab% nv, tab% lgRho, lgPr, lgRho, dlgRho, ierr)
		if (ierr /= 0) then
			write(error_unit,'(a,i3)') routine_name//': ierr = ',ierr
		end if
		call interp_value_and_slope(tab% lgP, tab% nv, tab% lgEps, lgPr, lgEps, dlgEps, ierr)
		if (ierr /= 0) then
			write(error_unit,'(a,i3)') routine_name//': ierr = ',ierr
		end if
		call interp_value_and_slope(tab% lgP, tab% nv, tab% Xprot, lgPr, Xprot, dXprot, ierr)
		if (ierr /= 0) then
			write(error_unit,'(a,i3)') routine_name//': ierr = ',ierr
		end if
	end subroutine dStar_core_get_results
	
end module dStar_core_lib
