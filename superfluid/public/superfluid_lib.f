module superfluid_lib
    use superfluid_def
    
contains

    subroutine sf_startup(datadir,ierr)
        use exceptions_lib
        character(len=*), intent(in) :: datadir
        integer, intent(out) :: ierr
        type(alert) :: already_initialized = alert( &
        &   scope='sf_startup',message='module already initialized')

        if (sf_is_initialized) then
            call already_initialized% report
            return
        end if
        sf_datadir = trim(datadir)//'/Tc_data'
        sf_scale = 1.0_dp
        sf_is_initialized = .true.
        ierr = 0
    end subroutine sf_startup

    subroutine sf_shutdown()
        integer :: gap_id
        do gap_id = 1,max_number_sf_types
            call sf_free_one(gap_id)
        end do
        sf_is_initialized = .false.
    end subroutine sf_shutdown

    subroutine sf_load_gaps(p1pre,n1pre,n3pre,ierr)
        use exceptions_lib
        character(len=*), intent(in) :: p1pre,n1pre,n3pre
        integer, intent(out) :: ierr
        type(alert) :: uninitialized=alert(scope='sf_load_gaps', &
        &   message='module is unitialized')
        
        if (.not. sf_is_initialized) then
            ierr = -1
            call uninitialized% report
            return
        end if
        ierr = 0
        call sf_load_one(p1pre,proton_1S0,ierr)
        call sf_load_one(n1pre,neutron_1S0,ierr)
        call sf_load_one(n3pre,neutron_3P2,ierr)
    end subroutine sf_load_gaps

    subroutine sf_load_one(prefix,type,ierr)
        use load_sf, only : load_sf_table
        character(len=*), intent(in) :: prefix
        integer, intent(in) :: type
        integer, intent(out) :: ierr
        call load_sf_table(prefix,type,ierr)
    end subroutine sf_load_one

    subroutine sf_free_one(id)
        use load_sf, only : free_sf_table
        integer, intent(in) :: id
        type(sf_table_type), pointer :: tab
        tab => sf_tables(id)
        call free_sf_table(tab)
    end subroutine sf_free_one
    
    subroutine sf_get_results(kFp,kFn,Tc)
        real(dp), intent(in) :: kFp, kFn
        real(dp), dimension(max_number_sf_types), intent(out) :: Tc
        
        Tc = 0.0_dp
        if (sf_tables(proton_1S0)% is_loaded) Tc(proton_1S0) = sf_get_Tc(kFp, proton_1S0)
        if (sf_tables(neutron_1S0)% is_loaded) Tc(neutron_1S0) = sf_get_Tc(kFn, neutron_1S0)
        if (sf_tables(neutron_3P2)% is_loaded) Tc(neutron_3P2) = sf_get_Tc(kFn, neutron_3P2)
    end subroutine sf_get_results
    
    function sf_get_Tc(kF,id) result(Tc)
        use exceptions_lib
        use interp_1d_lib, only : interp_values
        real(dp), intent(in) :: kF
        integer, intent(in) :: id
        real(dp) :: Tc
        integer, parameter :: Np = 1
        real(dp), dimension(Np) :: k, T
        integer :: N, ierr
        type(sf_table_type), pointer :: tab
        type(failure) :: interpolation_error=failure(scope='sf_get_Tc', &
        &   message='interpolation error')

        tab => sf_tables(id)
        N = size(tab% kF)
        k = kF
        if ((kF < tab% kF_min) .or. (kF > tab% kF_max)) then
            Tc = 0.0
            return
        end if
        call interp_values(tab% kF, N, tab% f, Np, k, T, ierr)
        if (interpolation_error% raised(ierr)) then
            Tc = 0.0_dp
            return
        end if
            Tc = max(T(1),0.0_dp)*sf_scale(id)
    end function sf_get_Tc

end module superfluid_lib
