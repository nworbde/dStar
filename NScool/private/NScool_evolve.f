module NScool_evolve
    use NScool_def
    
contains
    subroutine get_nuclear_heating(s, ierr)
        use constants_def
        type(NScool_info), pointer :: s
        integer, intent(out) :: ierr
        real(dp) :: Ptop, Pbot, Q
        
        ierr = 0
        s% enuc = 0.0_dp
        
        ! outer crust
        Ptop = 10.0**s% lgP_min_heating_outer
        Pbot = 10.0**s% lgP_max_heating_outer
        Q = s% Q_heating_outer
        call do_one
        
        ! inner crust
        Ptop = 10.0**s% lgP_min_heating_inner
        Pbot = 10.0**s% lgP_max_heating_inner
        Q = s% Q_heating_inner
        call do_one
        
    contains
        subroutine do_one()
            real(dp) :: M_norm
            M_norm = sum(s% dm, s% P >= Ptop .and. s% P <= Pbot)
            where (s% P >= Ptop .and. s% P <= Pbot)
                s% enuc = s% Mdot * Q * mev_to_ergs * avogadro/M_norm
            end where
        end subroutine do_one
    end subroutine get_nuclear_heating    

end module NScool_evolve
