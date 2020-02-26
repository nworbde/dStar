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
        rhoC = s% eos_pasta_transition_in_fm3 * amu * density_n
        
        where(s% rho_bar > rhoC)
            s% ionic_bar% Q = s% extra_real_controls(2)
        end where
    end subroutine alt_Qimp
    
    subroutine alt_heating(id,ierr)
        use constants_def
        use NScool_def
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        type(NScool_info), pointer :: s
        real(dp) :: rhoC, rhoC1, Qalt, M_norm
        
        call get_NScool_info_ptr(id,s,ierr)
        if (ierr /= 0) return

        Qalt = s% extra_real_controls(1)
        rhoC = s% extra_real_controls(3)
        rhoC1= s% extra_real_controls(4)
        
        M_norm = sum(s% dm, s% rho >= rhoC .and. s% rho < rhoC1)
        where (s% rho >= rhoC .and. s% rho < rhoC1)
            s% enuc = s% Mdot * Qalt * mev_to_ergs * avogadro/M_norm
        end where

    end subroutine alt_heating

end module alt_micro
