module dStar_crust_mod
    use dStar_crust_def
    integer, parameter :: crust_filename_length=128
contains
    
    subroutine do_load_crust_table(prefix,eos_handle,Tref,ierr)
        use, intrinsic :: iso_fortran_env, only: error_unit
        character(len=*), intent(in) :: prefix
        integer, intent(in) :: eos_handle
        real(dp), intent(in) :: Tref
        integer, intent(out) :: ierr
        logical, parameter :: dbg = .FALSE.
        type(crust_table_type), pointer :: tab
        real(dp), pointer, dimension(:,:) :: lgRho_val    
        character(len=crust_filename_length) :: table_name, cache_filename
        logical :: have_cache
        integer :: unitno
        
        tab => crust_table
        ! if the table is already allocated, issue a warning and scrub the table
        if (tab% is_loaded) then
            write(error_unit,'(a)') 'do_load_crust_table: overwriting already loaded table'
            call do_free_crust_table(tab)
        end if

        call generate_crust_filename(prefix,Tref,table_name)
        cache_filename = trim(crust_datadir)//'/cache/'//trim(table_name)//'.bin'
        inquire(file=cache_filename,exist=have_cache)
        if (have_cache) then
            call do_read_crust_cache(cache_filename,tab,ierr)
            if (ierr == 0) return
        end if
        
        ! if we don't have the table, or could not load it, then generatue a 
        ! new one and write to cache
        tab% nv = crust_default_number_table_points
        tab% lgP_min = crust_default_lgPmin
        tab% lgP_max = crust_default_lgPmax
        call do_generate_crust_table(prefix,eos_handle,Tref,tab)
        tab% is_loaded = .TRUE.
        
        lgRho_val(1:4,1:tab% nv) => tab% lgRho(1:4*tab% nv)
        
        ! write informative message about range of table
        if (dbg) then
            write(error_unit,'(a,2f8.3)') 'do_load_crust_table: lgNb min, max = ', &
            &   10.0**(minval(lgRho_val(1,:)-log10(amu))-39.0), 10.0**(maxval(lgRho_val(1,:)-log10(amu))-39.0)
            write(error_unit,'(t21,a,2f8.3)') 'lgP min, max = ', &
            &   tab% lgP_min, tab% lgP_max
        end if
        
        if (.not.have_cache) then
            call do_write_crust_cache(cache_filename,tab,ierr)
        end if
        
    end subroutine do_load_crust_table
    
    subroutine do_generate_crust_table(prefix,eos_handle,Tref,tab)
        use constants_def, only: avogadro
        use nucchem_def, only: composition_info_type
        use hz90
        use interp_1d_def
        use interp_1d_lib
        character(len=*), intent(in) :: prefix
        integer, intent(in) :: eos_handle
        real(dp), intent(in) :: Tref
        type(crust_table_type), pointer :: tab
        real(dp), dimension(:), pointer :: work=>null()
        real(dp), pointer, dimension(:,:) :: lgRho_val, lgEps_val
        real(dp) :: lgPmin, delta_lgP
        integer :: N, i, ierr
        integer :: ncharged
        integer, dimension(HZ90_number) :: charged_ids
        real(dp), dimension(:), allocatable :: lgP, Xneut, lgRho, lgEps
        real(dp), dimension(:,:), allocatable :: Yion
        type(composition_info_type), dimension(:), allocatable :: ion_info
            
        lgPmin = tab% lgP_min
        delta_lgP = tab% lgP_max - lgPmin
        N = tab% nv
        allocate(lgP(N), Xneut(N), lgRho(N), lgEps(N), Yion(HZ90_number,N), ion_info(N))
        
        lgP = [ (lgPmin + real(i-1,dp)*(delta_lgP)/real(N-1,dp), i = 1, N)]
        
        call do_make_crust(lgP,Yion,Xneut,charged_ids,ncharged,ion_info)        
        call find_densities(eos_handle,lgP,lgRho,lgEps, &
            & Yion,ncharged,charged_ids,ion_info,Tref)
        
        call do_allocate_crust_table(tab, N, ierr)
        lgRho_val(1:4,1:N) => tab% lgRho(1:4*N)
        lgEps_val(1:4,1:N) => tab% lgEps(1:4*N)
        
        tab% lgP = lgP
        lgRho_val(1,:) = lgRho
        lgEps_val(1,:) = lgEps
                
        allocate(work(tab% nv*pm_work_size))
        call interp_pm(tab% lgP, tab% nv, tab% lgRho, pm_work_size, work, &
        &   'do_generate_crust_table: Rho', ierr)
        call interp_pm(tab% lgP, tab% nv, tab% lgEps, pm_work_size, work, &
        &   'do_generate_crust_table: Eps', ierr)
        deallocate(work)
        deallocate(lgP, Xneut, lgRho, lgEps, Yion, ion_info)

    end subroutine do_generate_crust_table
    
    subroutine do_read_crust_cache(cache_filename,tab,ierr)
        character(len=*), intent(in) :: cache_filename
        type(crust_table_type), pointer :: tab
        integer, intent(out) :: ierr
        integer :: unitno, n
        
        open(newunit=unitno,file=trim(cache_filename), &
        &   action='read',status='old',form='unformatted', iostat=ierr)
        if (failure('opening'//trim(cache_filename),ierr)) return
        
        read(unitno) n
        call do_allocate_crust_table(tab,n,ierr)
        if (failure('allocating table',ierr)) then
            close(unitno)
            return
        end if
        
        read(unitno) tab% lgP(:)
        read(unitno) tab% lgRho
        read(unitno) tab% lgEps
        close(unitno)
        
        tab% lgP_min = minval(tab% lgP)
        tab% lgP_max = maxval(tab% lgP)
        tab% is_loaded = .TRUE.
    end subroutine do_read_crust_cache
    
    subroutine do_write_crust_cache(cache_filename,tab,ierr)
        character(len=*), intent(in) :: cache_filename
        type(crust_table_type), pointer :: tab
        integer, intent(out) :: ierr
        integer :: unitno
        
        ierr = 0
        if (tab% nv == 0) return
        open(newunit=unitno, file=trim(cache_filename),action='write', &
        &   form='unformatted',iostat=ierr)
        if (failure('opening '//trim(cache_filename),ierr)) return
        
        write(unitno) tab% nv
        write(unitno) tab% lgP
        write(unitno) tab% lgRho
        write(unitno) tab% lgEps
        close(unitno)
    end subroutine do_write_crust_cache
    
    subroutine do_allocate_crust_table(tab,n,ierr)
        type(crust_table_type),pointer :: tab
        integer, intent(in) :: n
        integer, intent(out) :: ierr
                
        allocate(tab% lgP(n), tab% lgRho(4*n), tab% lgEps(4*n), stat=ierr)
        if (ierr /= 0) return
        tab% nv = n
    end subroutine do_allocate_crust_table
    
    subroutine generate_crust_filename(prefix,Tref,filename)
        ! naming convention for flies is prefix_ttt
        ! where ttt = 100*log10(Tref), to 3 significant digits
        character(len=*), intent(in) :: prefix
        real(dp), intent(in) :: Tref
        character(len=crust_filename_length), intent(out) :: filename
        
        write (filename,'(a,"_",i0.3)') trim(prefix), int(100.0*log10(Tref))
    end subroutine generate_crust_filename

    subroutine do_free_crust_table(tab)
        type(crust_table_type), pointer :: tab
        tab% nv = 0
        tab% lgP_min = 0.0
        tab% lgP_max = 0.0
        if (allocated(tab% lgP)) deallocate(tab% lgP)
        if (associated(tab% lgRho)) deallocate(tab% lgRho)
        if (associated(tab% lgEps)) deallocate(tab% lgEps)
        nullify(tab% lgRho)
        nullify(tab% lgEps)
        tab% is_loaded = .FALSE.
    end subroutine do_free_crust_table
        
    function failure(msg,ierr)
        use, intrinsic :: iso_fortran_env, only: error_unit
        character(len=*), intent(in) :: msg
        integer, intent(in) :: ierr
        logical :: failure
        
        failure = (ierr /= 0)
        if (failure) then
            write(error_unit,*) trim(msg),': ierr = ',ierr
        end if
    end function failure
end module dStar_crust_mod
