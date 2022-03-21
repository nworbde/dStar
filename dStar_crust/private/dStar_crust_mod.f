module dStar_crust_mod
    use dStar_crust_def
    integer, parameter :: crust_filename_length=128

contains
    
    subroutine do_load_crust_table(prefix,eos_handle,Tref,ierr)
        use math_lib
        use exceptions_lib
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
        type(alert) :: overwrite=alert(scope='do_load_crust_table', &
        &   message='overwriting previously loaded table')
        character(len=128) :: alert_msg
        type(alert) :: status=alert(scope='do_load_crust_table')
        type(warning) :: cache_write_failure=warning(scope='do_load_crust_table', &
        &   message='unable to write cache')
        type(failure) :: crust_table_failure=failure(scope='do_load_crust_table', &
            message='unable to generate crust table')
        
        ierr = 0
        tab => crust_table
        ! if the table is already allocated, issue a warning and scrub the table
        if (tab% is_loaded) then
            call overwrite% report
            call do_free_crust_table(tab)
        end if

        ! attempt to read from a cache file
        call generate_crust_filename(prefix,Tref,table_name)
        cache_filename = trim(crust_datadir)//'/cache/'//trim(table_name)//'.bin'
        inquire(file=cache_filename,exist=have_cache)
        if (have_cache) then
            call status% report('loading crust table '//trim(cache_filename))
            call do_read_crust_cache(cache_filename,tab,ierr)
            if (ierr == 0) return
        end if
        
        ! if we don't have the table, or could not load it, then generatue a 
        ! new one and write to cache
        ierr = 0
        write(alert_msg,'(a,es9.2,a)') 'generating crust table using composition profile '//trim(prefix)// &
        &   ' at T = ',Tref,' K'
        call status% report(alert_msg)
        call do_generate_crust_table(prefix,eos_handle,Tref, tab,ierr)
        if (crust_table_failure% raised(ierr)) return
        tab% is_loaded = .TRUE.

        ! write informative message about range of table
        if (dbg) then
            write(alert_msg,'(a,2f8.3)') 'lgNb min, max = ', &
            &   exp10(minval(lgRho_val(1,:)-log10(amu))-39.0_dp), &
            &   exp10(maxval(lgRho_val(1,:)-log10(amu))-39.0_dp)
            call status% report(alert_msg)
            write(alert_msg,'(t21,a,2f8.3)') 'lgP min, max = ', &
            &   tab% lgP_min, tab% lgP_max
            call status% report(alert_msg)
        end if

        if (.not.have_cache) then
            call do_write_crust_cache(cache_filename,tab,ierr)
            if (cache_write_failure% raised(ierr)) return
        end if
    end subroutine do_load_crust_table
    
    subroutine do_generate_crust_table(prefix,eos_handle,T,tab,ierr)
        use exceptions_lib
        use constants_def, only: avogadro
        use nucchem_def
        use nucchem_lib
        use interp_1d_def
        use interp_1d_lib
        use composition_handler
        character(len=*), intent(in) :: prefix
        integer, intent(in) :: eos_handle
        real(dp), intent(in) :: T
        type(crust_table_type), pointer :: tab
        integer, intent(out) :: ierr
        ! for the composition table
        character(len=crust_filename_length) :: composition_filename
        integer :: N, Nisos
        real(dp), dimension(:), allocatable :: lgP
        character(len=iso_name_length), dimension(:), allocatable :: isos
        real(dp), dimension(:,:), allocatable :: Y
        ! interpolation storage
        real(dp), dimension(:), pointer :: work=>null()
        real(dp), pointer, dimension(:,:) :: Y_val
        real(dp), pointer, dimension(:) :: Yptr
        real(dp), pointer, dimension(:,:) :: lgRho_val, lgEps_val
        real(dp), dimension(:), allocatable :: lgRho,lgEps
        real(dp), dimension(:,:), allocatable :: Yion
        type(composition_info_type), dimension(:), allocatable :: ion_info
        integer, dimension(:), allocatable :: charged_ids
        integer, dimension(:), allocatable :: indcs
        integer :: ncharged, i
        real(dp) :: Xsum
        character(len=*), parameter :: routine_name='do_generate_crust_table'
        type(failure) :: read_composition_error=failure(scope=routine_name, &
        &   message='unable to read composition')
        type(failure) :: allocation_error=failure(scope=routine_name)
        type(failure) :: interpolation_error=failure(scope=routine_name)
                
        ierr = 0
        composition_filename = trim(crust_datadir)//'/'//trim(prefix)//'.bin'
        call read_composition_cache(composition_filename,N,Nisos,isos,lgP,Y,ierr)
        if (read_composition_error% raised(ierr)) return

        tab% nv = N
        tab% nisos = Nisos
        tab% T = T
        tab% lgP_min = minval(lgP)
        tab% lgP_max = maxval(lgP)
        allocate(lgRho(N), lgEps(N), Yion(Nisos,N), ion_info(N), indcs(Nisos), charged_ids(Nisos), stat=ierr)
        if (allocation_error% raised(ierr,'unable to allocate storage')) return
        
        ! set up markers for the nuclei
        indcs = [(get_nuclide_index(adjustl(isos(i))),i=1,Nisos)]

        ! set the composition moments and densities
        do i = 1, N
            call compute_composition_moments(Nisos, indcs, Y(:,i), &
            &   ion_info(i), Xsum, ncharged, charged_ids, Yion(:,i), &
            &   abunds_are_mass_fractions=.FALSE., exclude_neutrons=.TRUE.)
        end do
        call find_densities(eos_handle,lgP,lgRho,lgEps,Yion,ncharged,charged_ids,ion_info,tab% T)
        
        call do_allocate_crust_table(tab, N, Nisos, ierr)
        if (allocation_error% raised(ierr,'unable to allocate crust table')) return
        
        ! copy network information to table
        tab% network = isos

        ! construct interpolants
        tab% lgP = lgP
        allocate(work(tab% nv*pm_work_size))
        
        ! composition
        do i = 1, Nisos
            Y_val(1:4,1:N) => tab% Y(1:4*N,i)
            Y_val(1,:) = Y(i,:)
            Yptr(1:4*N) => tab% Y(1:4*N,i)
            call interp_pm(tab% lgP, tab% nv, Yptr, pm_work_size, work, &
            &   'do_generate_default_crust_table: Y',ierr)
            if (interpolation_error% raised(ierr,trim(nuclib% name(indcs(i))))) return
        end do

        ! eos
        lgRho_val(1:4,1:N) => tab% lgRho(1:4*N)
        lgEps_val(1:4,1:N) => tab% lgEps(1:4*N)
        lgRho_val(1,:) = lgRho
        lgEps_val(1,:) = lgEps

        call interp_pm(tab% lgP, tab% nv, tab% lgRho, pm_work_size, work, &
        &   routine_name//': Rho', ierr)
        if (interpolation_error% raised(ierr,'unable to interpolate lgRho')) return
        call interp_pm(tab% lgP, tab% nv, tab% lgEps, pm_work_size, work, &
        &   routine_name//': Eps', ierr)
        if (interpolation_error% raised(ierr,'unable to interpolate lgEps')) return
        
        deallocate(work)
        deallocate(lgP, lgRho, lgEps, Yion, ion_info, indcs, charged_ids, isos, Y)

    end subroutine do_generate_crust_table
    
    subroutine do_read_crust_cache(cache_filename,tab,ierr)
        use exceptions_lib
        character(len=*), intent(in) :: cache_filename
        type(crust_table_type), pointer :: tab
        integer, intent(out) :: ierr
        integer :: unitno, nv, nisos
        type(failure) :: cache_error=failure(scope='do_read_crust_cache')
        
        open(newunit=unitno,file=trim(cache_filename), &
        &   action='read',status='old',form='unformatted', iostat=ierr)
        if (cache_error% raised(ierr,'opening'//trim(cache_filename))) &
        &   return
        
        read(unitno) nv
        read(unitno) nisos
        call do_allocate_crust_table(tab,nv,nisos,ierr)
        if (cache_error% raised(ierr,'allocating table')) then
            close(unitno)
            return
        end if
        
        tab% nv = nv
        tab% nisos = nisos
        read(unitno) tab% network
        read(unitno) tab% T
        read(unitno) tab% lgP
        read(unitno) tab% Y
        read(unitno) tab% lgRho
        read(unitno) tab% lgEps
        close(unitno)
        
        tab% lgP_min = minval(tab% lgP)
        tab% lgP_max = maxval(tab% lgP)
        tab% is_loaded = .TRUE.
    end subroutine do_read_crust_cache
    
    subroutine do_write_crust_cache(cache_filename,tab,ierr)
        use exceptions_lib
        character(len=*), intent(in) :: cache_filename
        type(crust_table_type), pointer :: tab
        integer, intent(out) :: ierr
        integer :: unitno
        type(failure) :: cache_error=failure(scope='do_write_crust_cache')
        
        ierr = 0
        if (tab% nv == 0 .or. tab% nisos == 0) return
        open(newunit=unitno, file=trim(cache_filename),action='write', &
        &   form='unformatted',iostat=ierr)
        if (cache_error% raised(ierr,'opening '//trim(cache_filename))) &
        &   return
        
        write(unitno) tab% nv
        write(unitno) tab% nisos
        write(unitno) tab% network
        write(unitno) tab% T
        write(unitno) tab% lgP
        write(unitno) tab% Y
        write(unitno) tab% lgRho
        write(unitno) tab% lgEps
        close(unitno)
    end subroutine do_write_crust_cache
    
    subroutine do_allocate_crust_table(tab,n,nisos,ierr)
        use exceptions_lib
        type(crust_table_type), pointer :: tab
        integer, intent(in) :: n,nisos
        integer, intent(out) :: ierr
        type(failure) :: allocation_error=failure(scope='do_allocate_crust_table', &
        &   message='unable to allocate crust table')

        tab% nv = 0
        tab% nisos = 0
        allocate(tab% lgP(n),tab% lgEps(4*n), tab% lgRho(4*n), &
        &   tab% network(nisos), tab% Y(4*n,nisos), stat=ierr)
        if (allocation_error% raised(ierr)) return                
        tab% nv = n
        tab% nisos = nisos
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
        tab% nisos = 0
        tab% lgP_min = 0.0
        tab% lgP_max = 0.0
        tab% T = 0.0
        
        if (allocated(tab% lgP)) deallocate(tab% lgP)
        if (allocated(tab% network)) deallocate(tab% network)
        if (associated(tab% Y)) deallocate(tab% Y)
        if (associated(tab% lgRho)) deallocate(tab% lgRho)
        if (associated(tab% lgEps)) deallocate(tab% lgEps)
        nullify(tab% Y)
        nullify(tab% lgRho)
        nullify(tab% lgEps)

        tab% is_loaded = .FALSE.
    end subroutine do_free_crust_table

    subroutine find_densities(eos_handle,lgP,lgRho,lgEps, &
            & Yion,ncharged,charged_ids,ionic,Tref)
        use math_lib
        use exceptions_lib
        use constants_def
        use nucchem_def
        use num_lib

        integer, intent(in) :: eos_handle
        real(dp), dimension(:), intent(in) :: lgP
        real(dp), dimension(:), intent(out) :: lgRho
        real(dp), dimension(:), intent(out) :: lgEps
        real(dp), dimension(:,:), intent(in) :: Yion
        integer, intent(in) :: ncharged
        integer, dimension(:), intent(in) :: charged_ids
        type(composition_info_type), dimension(:), intent(in) :: ionic
        real(dp), intent(in) :: Tref
        real(dp) :: Pfac
        real(dp), dimension(:), pointer :: rpar=>null()
        integer, dimension(:), pointer :: ipar=>null()
        integer :: lipar, lrpar
        integer :: i,Ntab
        real(dp) :: x1, x3, y1, y3, epsx, epsy, lgRho_guess
        integer :: imax, ierr
        type(warning) :: rootfind_error=warning(scope='find_densities')
        type(alert) :: rootfind_status=alert(scope='find_densities')
        character(len=128) :: error_message
        
        Pfac = 0.25_dp*(threepisquare)**onethird *hbar*clight*avo**(4.0*onethird)
        Ntab = size(lgP)
        imax = 20
        epsx = 1.0d-12
        epsy = 1.0d-12
        
        ! decide the size of the parameter arrays
        lipar = 2 + ncharged
        allocate(ipar(lipar))
        ipar(1) = eos_handle
        ipar(2) = ncharged
        ipar(3:ncharged+2) = charged_ids(1:ncharged)
        
        lrpar = ncharged + 11 + 4
        allocate(rpar(lrpar))

        ! last value of rpar is a guess for the density; if 0, will be calculated for relativistic electron gas
        rpar(lrpar) = 0.0
        do i = 1, Ntab
            ! stuff the composition information into the parameter array
            rpar(1:ncharged) = Yion(1:ncharged,i)
            rpar(ncharged+1) = ionic(i)% A
            rpar(ncharged+2) = ionic(i)% Z
            rpar(ncharged+3) = ionic(i)% Z53
            rpar(ncharged+4) = ionic(i)% Z2
            rpar(ncharged+5) = ionic(i)% Z73
            rpar(ncharged+6) = ionic(i)% Z52
            rpar(ncharged+7) = ionic(i)% ZZ1_32
            rpar(ncharged+8) = ionic(i)% Z2XoA2
            rpar(ncharged+9) = ionic(i)% Ye
            rpar(ncharged+10) = ionic(i)% Yn
            rpar(ncharged+11) = ionic(i)% Q
            rpar(ncharged+12) = lgP(i)
            rpar(ncharged+13) = Tref
            
            if (i > 1 .and. rpar(lrpar) /= 0.0) then
                lgRho_guess = lgRho(i-1) + (lgP(i)-lgP(i-1))/rpar(lrpar)
            else
                lgRho_guess = log10((10.0**lgP(i)/Pfac)**0.75 /ionic(i)% Ye)
            end if
            
            call look_for_brackets(lgRho_guess,0.05*lgRho_guess,x1,x3,match_density, &
            &   y1,y3,imax,lrpar,rpar,lipar,ipar,ierr)
            if (rootfind_error% raised(ierr)) then
                write(error_message,'(a,5f0.3)') &
                &   'unable to bracket root: ',lgP(i), x1, x3, y1, y3
                call rootfind_status% report(error_message)
                cycle
            end if
            
            lgRho(i) = safe_root_with_initial_guess(match_density,lgRho_guess,x1,x3,y1,y3, &
            &   imax,epsx,epsy,lrpar,rpar,lipar,ipar,ierr)
            if (rootfind_error% raised(ierr)) then
                write(error_message,'(a,5f0.3)') &
                &   'unable to converge: ',lgP(i), x1, x3, y1, y3
                call rootfind_status% report(error_message)
                cycle
            end if
            lgEps(i) = rpar(ncharged+14)
        end do
    end subroutine find_densities
    
    real(dp) function match_density(lgRho, dfdlgRho, lrpar, rpar, lipar, ipar, ierr)
        ! returns with ierr = 0 if was able to evaluate f and df/dx at x
        ! if df/dx not available, it is okay to set it to 0
        use ieee_arithmetic
        use exceptions_lib
        use constants_def
        use superfluid_def, only: max_number_sf_types
        use superfluid_lib, only: sf_get_results
        use nucchem_def
        use dStar_eos_def
        use dStar_eos_lib
        
        integer, intent(in) :: lrpar, lipar
        real(dp), intent(in) :: lgRho
        real(dp), intent(out) :: dfdlgRho
        integer, intent(inout), pointer :: ipar(:) ! (lipar)
        real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
        integer, intent(out) :: ierr
        integer :: eos_handle, ncharged
        type(composition_info_type) :: ionic
        integer, dimension(:), allocatable :: charged_ids
        real(dp), dimension(:), allocatable :: Yion
        real(dp), dimension(num_dStar_eos_results) :: res
        integer :: phase
        real(dp) :: chi, lgPwant, lgP, kFn, kFp, Tcs(max_number_sf_types)
        real(dp) :: rho, T, Eint
        type(assertion) :: lgP_is_num = assertion(scope='match_density', &
        &    message='pressure is a number')
        type(crust_eos_component), dimension(num_crust_eos_components) :: eos_components
       
        ierr = 0
        eos_handle = ipar(1)
        ncharged = ipar(2)
        allocate(charged_ids(ncharged),Yion(ncharged))
        charged_ids(:) = ipar(3:ncharged+2)
        
        Yion = rpar(1:ncharged)
        ionic% A = rpar(ncharged+1)
        ionic% Z = rpar(ncharged+2)
        ionic% Z53 = rpar(ncharged+3)
        ionic% Z2 = rpar(ncharged+4)
        ionic% Z73 = rpar(ncharged+5)
        ionic% Z52 = rpar(ncharged+6)
        ionic% ZZ1_32 = rpar(ncharged+7)
        ionic% Z2XoA2 = rpar(ncharged+8)
        ionic% Ye = rpar(ncharged+9)
        ionic% Yn = rpar(ncharged+10)
        ionic% Q = rpar(ncharged+11)
        
        rho = 10.0**lgRho
        T = rpar(ncharged+13)
        ! when searching for a root, the density can be > nuclear, so that chi > 1.0,
        ! which results in a bad neutron number density.
        ! we cap chi at the value appropriate a density of 1.0e14 g/cc.
        chi = nuclear_volume_fraction(min(rho,1.0e14_dp),ionic,default_nuclear_radius)

        kFp = 0.0_dp
        kFn = neutron_wavenumber(rho,ionic,chi)
        call sf_get_results(kFp,kFn,Tcs)
        call eval_crust_eos(eos_handle,rho,T,ionic,ncharged, &
        &    charged_ids,Yion,Tcs,res,phase,chi,eos_components)
        Eint = exp(res(i_lnE))
        lgPwant = rpar(ncharged+12)
        lgP = res(i_lnP)/ln10
        call lgP_is_num% assert(.not. ieee_is_nan(lgP))
           
        rpar(ncharged+14) = log10(rho* &
        &   (1.0+dot_product(Yion(:),nuclib%mass_excess(charged_ids))/amu_n + Eint/clight2))
        rpar(lrpar) = res(i_chiRho)
        match_density = lgP - lgPwant
        deallocate(charged_ids,Yion)
    end function match_density

end module dStar_crust_mod
