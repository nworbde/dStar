module dStar_crust_lib
    use dStar_crust_def
    
contains
    subroutine dStar_crust_startup(datadir, ierr)
        use exceptions_lib
        character(len=*), intent(in) :: datadir
        integer, intent(out) :: ierr
        type(alert) :: initialization=alert(scope='dStar_crust_startup', &
        &   message='already initialized')

        if (crust_is_initialized) then
            ierr = 1
            call initialization% report
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
    
    subroutine dStar_crust_load_table(prefix,eos_handle,Tref,ierr)
        use exceptions_lib
        use dStar_crust_mod, only : do_load_crust_table
        character(len=*), intent(in) :: prefix
        integer, intent(in) :: eos_handle
        real(dp), intent(in) :: Tref
        integer, intent(out) :: ierr
        type(failure) :: load_error=failure(scope='dStar_crust_load_table')

        call do_load_crust_table(prefix, eos_handle, Tref, ierr)
        if (load_error% raised(ierr)) return
    end subroutine dStar_crust_load_table
    
    subroutine dStar_crust_get_results(lgP,lgRho,dlgRho,lgEps,dlgEps,ierr)
        use exceptions_lib
        use nucchem_def
        use interp_1d_lib, only : interp_value_and_slope
        real(dp), intent(in) :: lgP
        real(dp), intent(out) :: lgRho,lgEps,dlgRho,dlgEps
        integer, intent(out) :: ierr
        real(dp) :: lgP_c
        type(crust_table_type), pointer :: tab
        character(len=*), parameter :: routine_name = 'dStar_crust_get_results'
        type(assertion) :: is_loaded=assertion(scope=routine_name, &
        &   message='table is not loaded')
        type(failure) :: interpolation_error=failure(scope=routine_name)
        
        tab => crust_table
        call is_loaded% assert(tab% is_loaded)

        ! clip lgP to table
        lgP_c = max(lgP,tab% lgP_min)
        lgP_c = min(lgP_c,tab% lgP_max)

        call interp_value_and_slope(tab% lgP, tab% nv, tab% lgRho, lgP_c, lgRho, dlgRho, ierr)
        if (interpolation_error% raised(ierr,message='unable to interpolate lgRho')) &
        &   return
        call interp_value_and_slope(tab% lgP, tab% nv, tab% lgEps, lgP_c, lgEps, dlgEps, ierr)
        if (interpolation_error% raised(ierr,message='unable to interpolate lgEps')) &
        &   return
    end subroutine dStar_crust_get_results
        
    subroutine dStar_crust_get_composition(lgP,Yion,ierr)
        use exceptions_lib
        use nucchem_def
        use interp_1d_lib, only : interp_value
        real(dp), intent(in) :: lgP
        real(dp), dimension(:), intent(out) :: Yion
        integer, intent(out) :: ierr
        real(dp) :: lgP_c
        integer :: i
        real(dp), dimension(:), pointer :: Yptr
        type(crust_table_type), pointer :: tab
        character(len=*), parameter :: routine_name = 'dStar_crust_get_composition'
        type(assertion) :: is_loaded=assertion(scope=routine_name, &
        &   message='table is not loaded')
        type(failure) :: interpolation_error=failure(scope=routine_name, &
        &   message='unable to interpolate composition')

        tab => crust_table
        call is_loaded% assert(tab% is_loaded)
        
        ! clip lgP to table
        lgP_c = max(lgP,tab% lgP_min)
        lgP_c = min(lgP_c,tab% lgP_max)
        
        do i = 1, tab% Nisos
            Yptr(1:4*tab% nv) => tab% Y(1:4*tab% nv,i)
            call interp_value(tab% lgP, tab% nv, Yptr, lgP_c, Yion(i), ierr)
            if (interpolation_error% raised(ierr)) return
        end do
    end subroutine dStar_crust_get_composition
    
    subroutine dStar_crust_get_composition_info(lgP,ncharged,charged_ids,Yion,Xneut,ion_info,ierr)
        use exceptions_lib
        use nucchem_def
        use nucchem_lib    
        real(dp), dimension(:), intent(in) :: lgP   ! length = nz
        integer, intent(out) :: ncharged
        integer, dimension(:), intent(out) :: charged_ids   ! (tab% Nisos)
        real(dp), dimension(:,:), intent(out) :: Yion   ! (tab% Nisos,nz)
        real(dp), dimension(:), intent(out) :: Xneut    ! (nz)
        type(composition_info_type), dimension(:), intent(out) :: ion_info  ! (nz)
        integer, intent(out) :: ierr
        type(crust_table_type), pointer :: tab
        real(dp), dimension(:), allocatable :: Y
        integer :: N, Nisos, i
        integer, dimension(:), allocatable :: indcs
        real(dp) :: Xsum
        type(failure) :: get_composition_error=failure(scope='dStar_crust_get_composition_info', &
        &   message='unable to find composition')
        
        tab => crust_table
        N = size(lgP)
        Nisos = tab% Nisos
        
        allocate(Y(Nisos),indcs(Nisos))
        indcs = [(get_nuclide_index(adjustl(tab% network(i))),i=1,Nisos)]
        
        do i = 1, N
            call dStar_crust_get_composition(lgP(i),Y(:),ierr)
            if (get_composition_error% raised(ierr)) return
            ! ensure abundances are positive definite and renormalize
            where (Y(:) < 0.0_dp) Y(:) = 0.0_dp
            call compute_composition_moments(Nisos, indcs, Y(:), &
            &   ion_info(i), Xsum, ncharged, charged_ids, Yion(:,i), &
            &   abunds_are_mass_fractions=.FALSE., exclude_neutrons=.TRUE., &
            &   renormalize_mass_fractions=.TRUE.)
        end do
        Xneut = ion_info% Yn
        deallocate(Y,indcs)
    end subroutine dStar_crust_get_composition_info
    
    function dStar_crust_get_composition_size()
        integer :: dStar_crust_get_composition_size
        type(crust_table_type), pointer :: tab
        tab => crust_table
        dStar_crust_get_composition_size = tab% Nisos
    end function dStar_crust_get_composition_size

end module dStar_crust_lib
