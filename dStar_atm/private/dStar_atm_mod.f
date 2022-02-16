module dStar_atm_mod
    use dStar_atm_def
    integer, parameter :: atm_filename_length=128
contains
    
    subroutine do_load_atm_table(prefix,grav,Plight,Pb,ierr)
        use exceptions_lib
        character(len=*), intent(in) :: prefix
        real(dp), intent(in) :: grav,Plight,Pb
        integer, intent(out) :: ierr
        type(atm_table_type), pointer :: tab
        character(len=atm_filename_length) :: table_name, cache_filename
        logical :: have_cache
        integer :: unitno
        type(alert) :: status=alert(scope='do_load_atm_table')
        
        tab => atm_table
        ! if the table is already allocated, issue a warning and scrub the table
        if (tab% is_loaded) then
            call status% report('overwriting already loaded table')
            call do_free_atm_table(tab)
        end if

        call generate_atm_filename(prefix,grav,Plight,Pb,table_name)
        cache_filename = trim(atm_datadir)//'/cache/'//trim(table_name)//'.bin'
        inquire(file=cache_filename,exist=have_cache)
        if (have_cache) then
            call do_read_atm_cache(cache_filename,tab,ierr)
            if (ierr == 0) return
        end if
        
        ! if we don't have the table, or could not load it, then generate a 
        ! new one and write to cache
        call status% report('generating table '//trim(table_name))
        tab% nv = atm_default_number_table_points
        tab% lgTb_min = atm_default_lgTbmin
        tab% lgTb_max = atm_default_lgTbmax
        call do_generate_atm_table(prefix,grav,Plight,Pb,tab)
        tab% is_loaded = .TRUE.
        
        if (.not.have_cache) then
            call do_write_atm_cache(cache_filename,tab,ierr)
        end if
    end subroutine do_load_atm_table
    
    subroutine do_generate_atm_table(prefix,grav,Plight,Pb,tab)
        use exceptions_lib
        use pcy97
        use bc09
        use interp_1d_def
        use interp_1d_lib
        character(len=*), intent(in) :: prefix
        real(dp), intent(in) :: grav,Plight,Pb
        type(atm_table_type), pointer :: tab
        real(dp), dimension(:), pointer :: work=>null()
        real(dp), pointer, dimension(:,:) :: lgTeff_val, lgflux_val
        real(dp) :: lgTbmin, delta_lgTb, lgTeffmin, delta_lgTeff
        integer :: N, i, ierr
        real(dp), dimension(:), allocatable :: lgTeff, lgflux, lgTb
        type(assertion) :: prefix_okay=assertion(scope='do_generate_atm_table', &
        &   message='atmosphere prefix is okay')
        type(assertion) :: atm_model_made=assertion(scope='do_generate_atm_table', &
        &   message='successfully generated atmosphere model')

        N = tab% nv
        allocate(lgTeff(N),lgflux(N),lgTb(N))
        select case (trim(prefix))
            
            case('pcy97')
            lgTbmin = tab% lgTb_min
            delta_lgTb = tab% lgTb_max - lgTbmin
            lgTb = [ (lgTbmin + real(i-1,dp)*(delta_lgTb)/real(N-1,dp),  &
            &   i = 1, N) ]
            call do_get_pcy97_Teff(grav, Plight, lgTb, lgTeff, lgflux)
            
            case('bc09')
            lgTeffmin = default_lgTeff_min
            delta_lgTeff = default_lgTeff_max - default_lgTeff_min
            lgTeff = [ (lgTeffmin + real(i-1,dp)*(delta_lgTeff)/real(N-1,dp), &
            &   i = 1, N) ]
            call do_get_bc09_Teff(grav, Plight, Pb, lgTb, lgTeff, lgflux, ierr)
            call atm_model_made% assert(ierr == 0)
            
            case default
            call prefix_okay% assert(.FALSE.)
        end select 

        call do_allocate_atm_table(tab, N, ierr)
        lgTeff_val(1:4,1:N) => tab% lgTeff(1:4*N)
        lgflux_val(1:4,1:N) => tab% lgflux(1:4*N)
        
        tab% lgTb = lgTb
        lgTeff_val(1,:) = lgTeff
        lgflux_val(1,:) = lgflux
        
        allocate(work(tab% nv*pm_work_size))
        call interp_pm(tab% lgTb, tab% nv, tab% lgTeff, pm_work_size, work, &
        &   'do_generate_atm_table: Teff', ierr)
        call interp_pm(tab% lgTb, tab% nv, tab% lgflux, pm_work_size, work, &
        &   'do_generate_atm_table: flux', ierr)
        deallocate(work)
        deallocate(lgTeff,lgflux,lgTb)

    end subroutine do_generate_atm_table
    
    subroutine do_read_atm_cache(cache_filename,tab,ierr)
        use exceptions_lib
        character(len=*), intent(in) :: cache_filename
        type(atm_table_type), pointer :: tab
        integer, intent(out) :: ierr
        integer :: unitno, n
        type(failure) :: opening_file_failure=failure(scope='do_read_atm_cache')
        type(failure) :: allocating_table_failure=failure(scope='do_read_atm_cache', &
        &   message='allocating atmosphere table')
        
        open(newunit=unitno,file=trim(cache_filename), &
        &   action='read',status='old',form='unformatted', iostat=ierr)
        if (opening_file_failure% raised(ierr,'opening'//trim(cache_filename))) return
        
        read(unitno) n
        call do_allocate_atm_table(tab,n,ierr)
        if (allocating_table_failure% raised(ierr)) then
            close(unitno)
            return
        end if
        
        read(unitno) tab% lgTb(:)
        read(unitno) tab% lgTeff
        read(unitno) tab% lgflux
        close(unitno)
        
        tab% lgTb_min = minval(tab% lgTb)
        tab% lgTb_max = maxval(tab% lgTb)
        tab% is_loaded = .TRUE.
    end subroutine do_read_atm_cache
    
    subroutine do_write_atm_cache(cache_filename,tab,ierr)
        use exceptions_lib
        character(len=*), intent(in) :: cache_filename
        type(atm_table_type), pointer :: tab
        integer, intent(out) :: ierr
        integer :: unitno
        type(assertion) :: table_exists=assertion(scope='do_write_atm_cache', &
        &   message='table exists')
        type(failure) :: opening_table_failure=failure(scope='do_write_atm_cache')
        
        ierr = 0
        call table_exists% assert(tab% nv > 0)
        
        open(newunit=unitno, file=trim(cache_filename),action='write', &
        &   form='unformatted',iostat=ierr)
        if (opening_table_failure% raised(ierr,'opening '//trim(cache_filename))) return
        
        write(unitno) tab% nv
        write(unitno) tab% lgTb
        write(unitno) tab% lgTeff
        write(unitno) tab% lgflux
        close(unitno)
    end subroutine do_write_atm_cache
    
    subroutine do_allocate_atm_table(tab,n,ierr)
        type(atm_table_type),pointer :: tab
        integer, intent(in) :: n
        integer, intent(out) :: ierr
        
        allocate(tab% lgTb(n),tab% lgTeff(4*n), tab% lgflux(4*n),stat=ierr)
        if (ierr /= 0) return
        tab% nv = n
    end subroutine do_allocate_atm_table
    
    subroutine generate_atm_filename(prefix,grav,Plight,Pb,filename)
        ! naming convention for flies is prefix_gggg_pppp_bbbb
        ! where gggg = 100*log10(g), pppp = 100*log10(Plight), 
        ! and bbbb = 100*log10(Pb), to 4 significant digits
        character(len=*), intent(in) :: prefix
        real(dp), intent(in) :: grav,Plight,Pb
        character(len=atm_filename_length), intent(out) :: filename
        
        write (filename,'(a,3("_",i0.4))') trim(prefix), &
                & int(100.0*log10(grav)), int(100.0*log10(Plight)), int(100.0*log10(Pb))
    end subroutine generate_atm_filename

    subroutine do_free_atm_table(tab)
        type(atm_table_type), pointer :: tab
        tab% nv = 0
        tab% lgTb_min = 0.0
        tab% lgTb_max = 0.0
        if (allocated(tab% lgTb)) deallocate(tab% lgTb)
        if (associated(tab% lgTeff)) deallocate(tab% lgTeff)
        if (associated(tab% lgflux)) deallocate(tab% lgflux)
        nullify(tab% lgTeff)
        nullify(tab% lgflux)
        tab% is_loaded = .FALSE.
    end subroutine do_free_atm_table
end module dStar_atm_mod
