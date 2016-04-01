module alt_micro

contains
    subroutine alt_Qimp(id,ierr)
        use constants_def
        use NScool_def
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        type(NScool_info), pointer :: s
        real(dp) :: rhoC
        
        call get_NScool_info_ptr(id,s,ierr)
        if (ierr /= 0) return
        rhoC = s% extra_real_controls(1)
        
        where(s% rho_bar > rhoC)
            s% ionic_bar% Q = s% extra_real_controls(2)
        end where
    end subroutine alt_Qimp
    
    subroutine alt_sf(id,kp,kn,Tc)
        use constants_def
        use superfluid_def
        use superfluid_lib
        use NScool_def
        
        integer, intent(in) :: id   ! for accessing parameters
        real(dp), intent(in) :: kp, kn  ! proton, neutron wavevectors (fm**-1)
        real(dp), dimension(max_number_sf_types), intent(out) :: Tc ! (K)
        type(NScool_info), pointer :: s
        real(dp) :: rhoC, kC, n
        integer :: ierr
        
        Tc = 1.0e10
        call get_NScool_info_ptr(id,s,ierr)
        if (ierr /= 0) return
        
        rhoC = s% extra_real_controls(1)
        n = rhoC*avogadro
        kC = (0.5*threepisquare*n)**onethird / cm_to_fm
        if (kn > kC) then
            Tc(neutron_1S0) = 0.0
        end if
    end subroutine alt_sf

end module alt_micro
