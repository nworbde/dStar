module dStar_core_mod
    use dStar_core_def
    integer, parameter :: core_filename_length=128
contains
    
    subroutine do_load_core_table(model,ierr)
        use iso_fortran_env, only: error_unit
        character(len=*), intent(in) :: model
        integer, intent(out) :: ierr
        type(core_table_type), pointer :: tab
        character(len=core_filename_length) :: table_name, cache_filename
        logical :: have_cache
        integer :: unitno
        
        tab => core_table
        ! if the table is already allocated, issue a warning and scrub the table
        if (tab% is_loaded) then
            write(error_unit,'(a)') 'do_load_core_table: overwriting already loaded table'
            call do_free_core_table(tab)
        end if

        call generate_core_filename(model,table_name)
        cache_filename = trim(core_datadir)//'/cache/'//trim(table_name)//'.bin'
        inquire(file=cache_filename,exist=have_cache)
        if (have_cache) then
            call do_read_core_cache(cache_filename,tab,ierr)
            if (ierr == 0) return
        end if
        
        ! if we don't have the table, or could not load it, then generate a 
        ! new one and write to cache
        write (error_unit, '(a)') 'generating table '//trim(table_name)//'...'
        tab% nv = core_default_number_table_points
        call do_generate_core_table(model,tab)
        tab% is_loaded = .TRUE.
        write (error_unit,'(a)') 'done'
        
        if (.not.have_cache) then
            call do_write_core_cache(cache_filename,tab,ierr)
        end if      
    end subroutine do_load_core_table
    
    subroutine do_generate_core_table(model,tab)
        use brown_skyrme
        use load_skyrme_parameters
        use interp_1d_def
        use interp_1d_lib
        character(len=*), intent(in) :: model
        type(core_table_type), pointer :: tab
        real(dp), dimension(:), pointer :: work=>null()
        real(dp), dimension(:,:), pointer :: lgRho_val, lgEps_val, Xprot_val
        real(dp) :: lgRhomin, delta_lgRho
        integer :: N, i, ierr
        real(dp), dimension(:), allocatable :: rho, eps, Xprot, P, muhat, mue, cs2
        
        N = tab% nv
        allocate(rho(N),eps(N),Xprot(N),P(N),muhat(N),mue(N),cs2(N))
        call load_skyrme_table(model,ierr)
        lgRhomin = log10(core_default_Rhomin)
        delta_lgRho = log10(core_default_Rhomax)-lgRhomin
        rho = [ (lgRhomin + real(i-1,dp)*delta_lgRho/real(N-1,dp), &
            & i=1,N) ]
        rho = 10.0_dp**rho
        call find_beta_eq(rho,Xprot,eps,P,muhat,mue,cs2,ierr)
        if (failure('fatal error: unable to generate core model',ierr)) stop

        call do_allocate_core_table(tab,N, ierr)
        lgRho_val(1:4,1:N) => tab% lgRho(1:4*N)
        lgEps_val(1:4,1:N) => tab% lgEps(1:4*N)
        Xprot_val(1:4,1:N) => tab% Xprot(1:4*N)
        
        tab% lgP = log10(P)
        lgRho_val(1,:) = log10(rho)
        lgEps_val(1,:) = log10(eps)
        Xprot_val(1,:) = Xprot
        
        allocate(work(tab% nv*pm_work_size))
        call interp_pm(tab% lgP, tab% nv, tab% lgRho, pm_work_size, work, &
        &   'do_generate_core_table: Rho', ierr)
        call interp_pm(tab% lgP, tab% nv, tab% lgEps, pm_work_size, work, &
        &   'do_generate_core_table: Eps', ierr)
        call interp_pm(tab% lgP, tab% nv, tab% Xprot, pm_work_size, work, &
        &   'do_generate_core_table: Xprot', ierr)
        deallocate(work)
        deallocate(rho,eps,Xprot,P,muhat,mue,cs2)
        tab% lgP_min = minval(tab% lgP)
        tab% lgP_max = maxval(tab% lgP)
    end subroutine do_generate_core_table
    
    subroutine do_read_core_cache(cache_filename,tab,ierr)
        character(len=*), intent(in) :: cache_filename
        type(core_table_type), pointer :: tab
        integer, intent(out) :: ierr
        integer :: unitno, n
        
        open(newunit=unitno,file=trim(cache_filename), &
        &   action='read',status='old',form='unformatted', iostat=ierr)
        if (failure('opening'//trim(cache_filename),ierr)) return
        
        read(unitno) n
        call do_allocate_core_table(tab,n,ierr)
        if (failure('allocating table',ierr)) then
            close(unitno)
            return
        end if
        
        read(unitno) tab% lgP
        read(unitno) tab% lgRho
        read(unitno) tab% lgEps
        read(unitno) tab% Xprot
        close(unitno)
        
        tab% lgP_min = minval(tab% lgP)
        tab% lgP_max = maxval(tab% lgP)
        tab% is_loaded = .TRUE.
    end subroutine do_read_core_cache
    
    subroutine do_write_core_cache(cache_filename,tab,ierr)
        character(len=*), intent(in) :: cache_filename
        type(core_table_type), pointer :: tab
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
        write(unitno) tab% Xprot
        close(unitno)
    end subroutine do_write_core_cache
 
    subroutine do_allocate_core_table(tab,n,ierr)
        type(core_table_type),pointer :: tab
        integer, intent(in) :: n
        integer, intent(out) :: ierr
        
        allocate(tab% lgP(n),tab% lgRho(4*n), tab% lgEps(4*n),tab% Xprot(4*n), &
            & stat=ierr)
        if (ierr /= 0) return
        tab% nv = n
    end subroutine do_allocate_core_table
  
    subroutine generate_core_filename(model,filename)
        ! naming convention for flies is model_eos
        character(len=*), intent(in) :: model
        character(len=core_filename_length), intent(out) :: filename
        
        write (filename,'(a,"_eos")') trim(model)
    end subroutine generate_core_filename

    subroutine do_free_core_table(tab)
        type(core_table_type), pointer :: tab
        tab% nv = 0
        tab% lgP_min = 0.0
        tab% lgP_max = 0.0
        if (allocated(tab% lgP)) deallocate(tab% lgP)
        if (associated(tab% lgRho)) deallocate(tab% lgRho)
        if (associated(tab% lgEps)) deallocate(tab% lgEps)
        if (associated(tab% Xprot)) deallocate(tab% Xprot)
        nullify(tab% lgRho)
        nullify(tab% lgEps)
        nullify(tab% Xprot)
        tab% is_loaded = .FALSE.
    end subroutine do_free_core_table
        
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

end module dStar_core_mod
