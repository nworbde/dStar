module dStar_atm_lib
    use dStar_atm_def
    
contains
    subroutine dStar_atm_startup(datadir, ierr)
        use exceptions_lib
        character(len=*), intent(in) :: datadir
        integer, intent(out) :: ierr
        type(alert) :: already_initialized=alert(level=1,scope='dStar_atm_startup', &
        &   message='module already initialized')
        if (atm_is_initialized) then
            ierr = 1
            call already_initialized% report
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
        use exceptions_lib
        use dStar_atm_mod, only : do_load_atm_table
        character(len=*), intent(in) :: prefix
        real(dp), intent(in) :: grav,Plight,Pb
        integer, intent(out) :: ierr
        type(failure) :: load_atm_failure=failure(scope='dStar_atm_load_table')
        type(alert) :: status=alert(scope='dStar_atm_load_table')
        
        call status% report('loading atmosphere model '//prefix)
        call do_load_atm_table(prefix, grav, Plight, Pb, ierr)
        if (load_atm_failure% raised(ierr)) return
    end subroutine dStar_atm_load_table
    
    subroutine dStar_atm_get_results(lgTb,lgTeff,dlgTeff,lgflux,dlgflux,ierr)
        use exceptions_lib
        use interp_1d_lib, only : interp_value_and_slope
        real(dp), intent(in) :: lgTb
        real(dp), intent(out) :: lgTeff,lgflux,dlgTeff,dlgflux
        integer, intent(out) :: ierr
        real(dp) :: lgT
        type(atm_table_type), pointer :: tab
        character(len=*), parameter :: routine_name = 'dStar_atm_get_results'
        type(assertion) :: table_loaded=assertion(scope='dStar_atm_get_results', &
        &   message='table is loaded')
        type(alert) :: status=alert(scope='dStar_atm_get_results',level=0)
        tab => atm_table
        call table_loaded% assert(tab% is_loaded)

        ! clip lgTb to table
        lgT = min(max(lgTb,tab% lgTb_min),tab% lgTb_max)
        call interp_value_and_slope(tab% lgTb, tab% nv, tab% lgTeff, lgT, lgTeff, dlgTeff, ierr)
        if (ierr /= 0) call status% report('unable to interpolate lgTeff')
        call interp_value_and_slope(tab% lgTb, tab% nv, tab% lgflux, lgT, lgflux, dlgflux, ierr)
        if (ierr /= 0) call status% report('unable to interpolate flux')
    end subroutine dStar_atm_get_results

end module dStar_atm_lib
