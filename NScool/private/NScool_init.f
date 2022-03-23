module init
    use math_lib
    use NScool_def, only: get_NScool_info_ptr
    use NScool_private_def
    use constants_def
    use constants_lib
    
contains
    subroutine do_NScool_init(my_dStar_dir, ierr)
        use exceptions_lib
        character(len=*), intent(in) :: my_dStar_dir
        integer, intent(out) :: ierr
        type(assertion) :: init_okay = assertion(scope='do_NScool_init', &
        & message='unable to initialize constants. Is DSTAR_DIR set or specified in the inlist?')
        ierr = 0
        
        call math_init()
        call constants_init(my_dStar_dir,'',ierr)
        call init_okay% assert(ierr==0)
        call NScool_private_def_init()
    end subroutine do_NScool_init

    subroutine do_NScool_shutdown()
        use dStar_atm_lib
        use dStar_crust_lib
        use dStar_eos_lib
        use nucchem_lib
        use superfluid_lib
        
        call dStar_atm_shutdown
        call dStar_crust_shutdown
        call dStar_eos_shutdown
        call nucchem_shutdown
        call sf_shutdown
    end subroutine do_NScool_shutdown

    function alloc_NScool_data(ierr)
        integer, intent(out) :: ierr
        integer :: alloc_NScool_data
        type(NScool_info), pointer :: s
            
        ierr = 0
        alloc_NScool_data = NScool_private_alloc(ierr)
        if (ierr /= 0) return
        
        call get_NScool_info_ptr(alloc_NScool_data, s, ierr)
        if (ierr /= 0) return
        
        !zero information
        s% Lsurf = 0.0
        s% dlnLsdlnT = 0.0 ! change surface luminosity with temperature out outer zone
        s% Teff = 0.0      ! surface effective temperature
        s% Lcore = 0.0     ! core luminosity
        s% Psurf = 0.0     ! surface pressure
        s% Pcore = 0.0     ! pressure at core boundary
        s% Mcore = 0.0     ! mass of core
        s% Rcore = 0.0     ! radius of core
        s% grav = 0.0      ! surface gravity
        s% tsec = 0.0      ! current value of time in seconds
        s% dt = 0.0        ! value of timestep just taken, in seconds
        s% Mdot = 0.0      ! accretion rate measured at infinity [g/s]
        s% model = 0       ! counter that is incremented after each successful step
      
        ! information about the composition
        s% nisos = 0  ! number of isotopes
        s% ncharged = 0
        
!         nullify(s% iso_ids)
        nullify(s% charged_ids)
      
        ! administrative
        s% eos_handle = -1
     
        s% base_profile_filename = ''
        s% history_filename = ''
        s% profile_manifest_filename = ''
        
        ! zonal information
        s% target_resolution_lnP = 0.0               ! target (d lnP) of a zone
        s% nz = -1
      
        nullify(s% dm)
!         nullify(s% X)
        nullify(s% Yion)
        nullify(s% Xneut)
        nullify(s% ionic)
        nullify(s% P)
        nullify(s% lnT)
        nullify(s% T)
        nullify(s% ePhi)
        nullify(s% e2Phi)
        nullify(s% eLambda)
        nullify(s% rho)
        nullify(s% lnCp)
        nullify(s% Cp)
        nullify(s% dlnCp_dlnT)
        nullify(s% lnGamma)
        nullify(s% Gamma)
        nullify(s% lnGamma)
        nullify(s% enu)
        nullify(s% lnenu)
        nullify(s% dlnenu_dlnT)
        nullify(s% enuc)

        ! facial information
        nullify(s% m)
        nullify(s% area)
        nullify(s% L)
        nullify(s% dm_bar)
!         nullify(s% X_bar)
        nullify(s% Yion_bar)
        nullify(s% Xneut_bar)
        nullify(s% ionic_bar)
        nullify(s% P_bar)
        nullify(s% lnT_bar)
        nullify(s% T_bar)
        nullify(s% ePhi_bar)
        nullify(s% e2Phi_bar)
        nullify(s% eLambda_bar)
        nullify(s% rho_bar)
        nullify(s% Kcond)
        nullify(s% lnK)
        nullify(s% dlnK_dlnT)

        ! storage for interpolation of tabulated ln(enu), ln(Cp), ln(Kcond)
        s% n_tab = 0                       ! number of table points
        nullify(s% tab_lnT) ! (n_tab) junction points for interpolation; same for all grid points
        nullify(s% tab_lnenu)  ! (4*n_tab, nz) coefficients for ln(enu)
        nullify(s% tab_lnCp)  ! (4*n_tab, nz) coefficients for ln(Cp)
        nullify(s% tab_lnGamma)
        nullify(s% tab_lnK)   ! (4*n_tab, nz) coefficients for ln(Kcond)
        
        ! storage for Teff(t), F(t), Qb(t)
        nullify(s% t_monitor)
        nullify(s% Teff_monitor)
        nullify(s% Qb_monitor)
    end function alloc_NScool_data
    
    subroutine dealloc_NScool_data(id, ierr)
        use storage, only: free_NScool_data
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        call free_NScool_data(id, ierr)        
    end subroutine dealloc_NScool_data
    
    subroutine model_builder(id, ierr)
        use create_model
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        call do_create_crust_model(id, ierr)
    end subroutine model_builder
    
end module init
