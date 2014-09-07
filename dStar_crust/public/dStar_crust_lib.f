module dStar_crust_lib
    use dStar_crust_def
    
contains
	subroutine dStar_crust_startup(datadir, ierr)
        use, intrinsic :: iso_fortran_env, only: error_unit
		character(len=*), intent(in) :: datadir
		integer, intent(out) :: ierr
		if (crust_is_initialized) then
            ierr = 1
			write(error_unit,*) 'dStar_crust_startup: ', &
			&    'package alreading initialized'
			return
		end if
		crust_datadir = trim(datadir)//'/crust_data'
		crust_is_initialized = .TRUE.
        ierr = 0
	end subroutine dStar_crust_startup

	subroutine dStar_crust_shutdown()
		use dStar_crust_mod, only : do_free_crust_table
		type(crust_table_type), pointer :: tab
		tab => crust_table
		call do_free_crust_table(tab)
		crust_is_initialized = .FALSE.
	end subroutine dStar_crust_shutdown
    
	subroutine dStar_crust_free_table()
		use dStar_crust_mod, only : do_free_crust_table
		type(crust_table_type), pointer :: tab
		tab => crust_table
		call do_free_crust_table(tab)
	end subroutine dStar_crust_free_table
    
	subroutine dStar_crust_load_table(prefix,Tref,ierr)
		use iso_fortran_env, only : error_unit
		use dStar_crust_mod, only : do_load_crust_table
		character(len=*), intent(in) :: prefix
		real(dp), intent(in) :: Tref
		integer, intent(out) :: ierr
		
		call do_load_crust_table(prefix, Tref, ierr)
		if (ierr /= 0) then
			write(error_unit,'(a,i3)') 'dStar_crust_load_one: ierr = ',ierr
		end if
	end subroutine dStar_crust_load_table
    
	subroutine dStar_crust_get_results(lgP,lgRho,dlgRho,lgNb,dlgNb,ierr)
		use iso_fortran_env, only : error_unit
		use interp_1d_lib, only : interp_value_and_slope
		real(dp), intent(in) :: lgP
		real(dp), intent(out) :: lgRho,lgNb,dlgRho,dlgNg
		integer, intent(out) :: ierr
		real(dp) :: lgP_c
		type(crust_table_type), pointer :: tab
        character(len=*), parameter :: routine_name = 'dStar_crust_get_results'
		tab => crust_table
		if (.not.tab% is_loaded) then
            ierr = -9
			write(error_unit,*) routine_name//': table is not loaded'
			return
		end if

		! clip lgP to table
		lgP_c = max(lgP,tab% lgP_min)
		lgT_c = min(lgP,tab% lgP_max)
		call interp_value_and_slope(tab% lgP, tab% nv, tab% lgRho, lgP_c, lgRho, dlgRho, ierr)
		if (ierr /= 0) then
			write (error_unit,'(a,i3)') routine_name//': ierr = ',ierr
		end if
		call interp_value_and_slope(tab% lgP, tab% nv, tab% lgNb, lgP_c, lgNb, dlgNb, ierr)
		if (ierr /= 0) then
			write (error_unit,'(a,i3)') routine_name//': ierr = ',ierr
		end if
	end subroutine dStar_crust_get_results

end module dStar_crust_lib
