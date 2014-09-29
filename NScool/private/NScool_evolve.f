module NScool_evolve
    use NScool_def
    
    integer, parameter :: i_id = 1
    integer, parameter :: num_deriv_ipar = 1, num_deriv_rpar = 0
    
    
contains
    
    subroutine get_derivatives(n, x, h, y, f, lrpar, rpar, lipar, ipar, ierr)
        use const_def, only: dp
        integer, intent(in) :: n, lrpar, lipar
        real(dp), intent(in) :: x, h
        real(dp), intent(inout) :: y(n) ! okay to edit y if necessary (e.g., replace negative values by zeros)
        real(dp), intent(out) :: f(n) ! dy/dx
        integer, intent(inout), pointer :: ipar(:) ! (lipar)
        real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
        integer, intent(out) :: ierr ! nonzero means retry with smaller timestep.
        type(NScool_info), pointer :: s
        
        ierr = 0
        call get_NScool_info_ptr(ipar(i_id), s, ierr)
        ! fatal erors
        if (ierr /= 0) then
            write (*,*) 'unable to acces NScool_info in get_derivatives'
            stop
        end if
        if (n /= s% nz) then
            ierr = -9
            write (*,*) 'wrong number of equations in solver'
            stop
        end if
        
        s% lnT(1:n) = y(1:n)
        s% T(1:n) = exp(s% lnT(1:n))
        call interpolate_temps(s)
        
        call get_coefficients(s, ierr)
        call evaluate_luminosity(s,ierr)

        f(1:n-1) = ((s% L(2:n)*s% e2Phi_bar(2:n) - s% L(1:n-1)*s% e2Phi_bar(1:n-1))/s% dm(1:n-1)/s% e2Phi(1:n-1)  &
        &   + s% enuc(1:n-1) - s% enu(1:n-1))*s% ePhi(1:n-1)/s% Cp(1:n-1)/s% T(1:n-1)
        f(n) = 0.0
        
    end subroutine get_derivatives
    
    subroutine evaluate_luminosity(s,ierr)
        use dStar_atm_lib
        type(NScool_info), pointer :: s
        integer, intent(out) :: ierr
        real(dp) :: lgTb, lgTeff, dlgTeff, lgflux, dlgflux
        
        lgTb = s% lnT_bar(1)/ln10
            
    	call dStar_atm_get_results(lgTb,lgTeff,dlgTeff,lgflux,dlgflux, ierr)
        if (ierr /= 0) return
        s% Lsurf = s% area(1) * 10.0**lgflux     ! emergent luminosity
        s% dlnLsdlnT = dlgflux
        s% Teff = 10.0**lgTeff     ! surface effective temperature
        s% L(1) = s% Lsurf        
        s% L(2:s% nz) = -s% area(2:s% nz)**2 *s% rho_bar(2:s% nz)*s% Kcond(2:s% nz)  &
        &   * (s% ePhi(1:s% nz-1)*s% T(1:s% nz-1) - s% ePhi(2:s% nz)*s% T(2:s% nz))  &
        &   / s% ePhi_bar(2:s% nz) / s% dm_bar(2:s% nz)
    end subroutine evaluate_luminosity
    
    subroutine interpolate_temps(s)
        type(NScool_info), pointer :: s
        
        s% T_bar(1) = s% T(1)
        s% T_bar(2:s% nz) = 0.5*(s% T(1:s% nz-1)*s% dm(2:s% nz) +  &
               & s% T(2:s% nz)*s% dm(1:s% nz-1))/s% dm_bar(2:s% nz)
        s% lnT_bar = log(s% T_bar)
    end subroutine interpolate_temps
    
    subroutine get_coefficients(s, ierr)
        use interp_1d_lib
        type(NScool_info), pointer :: s
        integer, intent(out) :: ierr
        integer :: iz
        real(dp), dimension(:), pointer :: lnKcond_interp, lnCp_interp, lnEnu_interp
        
        do iz = 1, s% nz
            
            lnKcond_interp(1:4*s% n_tab) => s% tab_lnK(1:4*s% n_tab, iz)
            lnCp_interp(1:4*s% n_tab) => s% tab_lnCp(1:4*s% n_tab, iz)
            lnEnu_interp(1:4*s% n_tab) => s% tab_lnEnu(1:4*s% n_tab, iz)
        
    		call interp_value_and_slope(s% tab_lnT, s% n_tab, lnKcond_interp, s% lnT_bar(iz), s% lnK(iz), s% dlnK_dlnT(iz), ierr)
            if (failure('lnK', iz)) return

            call interp_value_and_slope(s% tab_lnT, s% n_tab, lnCp_interp, s% lnT(iz), s% lnCp(iz), s% dlnCp_dlnT(iz), ierr)
            if (failure('lnCp',iz)) return
            
            call interp_value_and_slope(s% tab_lnT, s% n_tab, lnEnu_interp, s% lnT(iz), s% lnenu(iz), s% dlnenu_dlnt(iz), ierr)
            if (failure('lnEnu',iz)) return
            
            s% Kcond = exp(s% lnK)
            s% Cp = exp(s% lnCp)
            s% enu = exp(s% lnenu)
            
            call get_nuclear_heating(s, ierr)
            
        end do
        
        contains
        function failure(str,izone)
           character(len=*), intent(in) :: str
           integer, intent(in) :: izone
           logical :: failure
         
           if (ierr == 0) then
              failure = .FALSE.
              return
           end if
           write (*,*) 'FAILURE interpolating ',str,' zone ',iz
           failure = .TRUE.
        end function failure
            
    end subroutine get_coefficients
    
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
