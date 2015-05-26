module dStar_atm_lib
    use dStar_atm_def
    
contains
	subroutine dStar_atm_startup(datadir, ierr)
        use, intrinsic :: iso_fortran_env, only: error_unit
		character(len=*), intent(in) :: datadir
		integer, intent(out) :: ierr
		if (atm_is_initialized) then
            ierr = 1
			write(error_unit,*) 'dStar_atm_startup: ', &
			&    'package alreading initialized'
			return
		end if
		atm_datadir = trim(datadir)//'/atm_data'
		atm_is_initialized = .TRUE.
        ierr = 0
	end subroutine dStar_atm_startup

	subroutine dStar_atm_shutdown()
		use dStar_atm_mod, only : do_free_atm_table
		type(atm_table_type), pointer :: tab
		tab => atm_table
		call do_free_atm_table(tab)
		atm_is_initialized = .FALSE.
	end subroutine dStar_atm_shutdown
    
	subroutine dStar_atm_free_table()
		use dStar_atm_mod, only : do_free_atm_table
		type(atm_table_type), pointer :: tab
		tab => atm_table
		call do_free_atm_table(tab)
	end subroutine dStar_atm_free_table
    
	subroutine dStar_atm_load_table(prefix,grav,Plight,Pb,ierr)
		use iso_fortran_env, only : error_unit
		use dStar_atm_mod, only : do_load_atm_table
		character(len=*), intent(in) :: prefix
		real(dp), intent(in) :: grav,Plight,Pb
		integer, intent(out) :: ierr
		
        write (error_unit,'(a)') 'loading atmosphere table '//prefix
		call do_load_atm_table(prefix, grav, Plight, Pb, ierr)
		if (ierr /= 0) then
			write(error_unit,'(a,i3)') 'dStar_atm_load_one: ierr = ',ierr
		end if
        write (error_unit,'(a)') 'done'
	end subroutine dStar_atm_load_table
    
	subroutine dStar_atm_get_results(lgTb,lgTeff,dlgTeff,lgflux,dlgflux,ierr)
		use iso_fortran_env, only : error_unit
		use interp_1d_lib, only : interp_value_and_slope
		real(dp), intent(in) :: lgTb
		real(dp), intent(out) :: lgTeff,lgflux,dlgTeff,dlgflux
		integer, intent(out) :: ierr
		real(dp) :: lgT
		type(atm_table_type), pointer :: tab
        character(len=*), parameter :: routine_name = 'dStar_atm_get_results'
		tab => atm_table
		if (.not.tab% is_loaded) then
            ierr = -9
			write(error_unit,*) routine_name//': table is not loaded'
			return
		end if

		! clip lgTb to table
		lgT = max(lgTb,tab% lgTb_min)
		lgT = min(lgTb,tab% lgTb_max)
		call interp_value_and_slope(tab% lgTb, tab% nv, tab% lgTeff, lgT, lgTeff, dlgTeff, ierr)
		if (ierr /= 0) then
			write (error_unit,'(a,i3)') routine_name//': ierr = ',ierr
		end if
		call interp_value_and_slope(tab% lgTb, tab% nv, tab% lgflux, lgT, lgflux, dlgflux, ierr)
		if (ierr /= 0) then
			write (error_unit,'(a,i3)') routine_name//': ierr = ',ierr
		end if
	end subroutine dStar_atm_get_results

end module dStar_atm_lib
