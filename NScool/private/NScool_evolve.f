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
    
    subroutine get_jacobian(n, x, h, y, f, dfdy, ldfy, lrpar, rpar, lipar, ipar, ierr)
        use const_def, only: dp
        integer, intent(in) :: n, ldfy, lrpar, lipar
        real(dp), intent(in) :: x, h
        real(dp), intent(inout) :: y(n)
        real(dp), intent(out) :: f(n) ! dy/dx
        real(dp), intent(out) :: dfdy(ldfy, n)
        ! dense: dfdy(i, j) = partial f(i) / partial y(j)
        ! banded: dfdy(i-j+mujac+1, j) = partial f(i) / partial y(j)
           ! uses rows 1 to mljac+mujac+1 of dfdy.
           ! The j-th column of the square matrix is stored in the j-th column of the
           ! array dfdy as follows:
           ! dfdy(mujac+1+i-j, j) = partial f(i) / partial y(j)
           ! for max(1, j-mujac)<=i<=min(N, j+mljac)
        integer, intent(inout), pointer :: ipar(:) ! (lipar)
        real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
        integer, intent(out) :: ierr ! nonzero means terminate integration
        type(NScool_info), pointer :: s
        ! work arrays
        real(dp), dimension(n) :: CTinv, CTdminv, ArKdm
        real(dp) :: wplus(1:n-1), wminus(2:n), dLpdlnT(1:n-1), dLdlnT(2:n)
        
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
        
        ! call to derivatives ensures that L, T, Tbar and their ln's are set, as well as the coefficients
        call get_derivatives(n, x, h, y, f, lrpar, rpar, lipar, ipar, ierr)
        
        CTinv(1:n) = s% ePhi(1:n)/(s% Cp(1:n)*s% T(1:n))
        CTdminv(1:n) = CTinv(1:n)/s% e2Phi(1:n)/ s% dm(1:n)
        ArKdm(1:n) = s% area(1:n)**2 * s% rho_bar(1:n) * s% Kcond(1:n)/s% ePhi_bar(1:n)/s% dm_bar(1:n)
        
        wplus(1:n-1) = s% dm(2:n)*s% T(1:n-1)/(s% dm(2:n)*s% T(1:n-1) + s% dm(1:n-1)*s% T(2:n))
        dLpdlnT(1:n-1) = s% L(2:n)*s% dlnK_dlnT(2:n)*wplus(1:n-1) - ArKdm(2:n)*s% ePhi(1:n-1)*s% T(1:n-1)
        
        wminus(2:n) = s% dm(1:n-1)*s% T(2:n)/(s% dm(2:n)*s% T(1:n-1)+ s% dm(1:n-1)*s% T(2:n))
        dLdlnT(2:n) = s% L(2:n)*s% dlnK_dlnT(2:n)*wminus(2:n) + Arkdm(2:n)*s% ePhi(2:n)*s% T(2:n)
        
        dfdy = 0.0_dp
        
        ! J(k-1,k) : 2 <= k <= n
        dfdy(1,2:n) = CTdminv(1:n-1)*s% e2Phi_bar(2:n)*dLdlnT(2:n)
        
        ! J(k,k) : 1 <= k <= n
        ! cell-based quantities
        dfdy(2,1:n) = -f(1:n)*(1.0_dp+s% dlnCp_dlnT(1:n)) - CTinv(1:n)*s% enu(1:n)*s% dlnenu_dlnT(1:n)
        ! interior points
        dfdy(2,2:n-1) = dfdy(2,2:n-1) + CTdminv(2:n-1)*(s% e2Phi_bar(3:n)*dLpdlnT(2:n-1) &
        &   - s% e2Phi_bar(2:n-1)*dLdlnT(2:n-1))
        ! surface point
        dfdy(2,1) = dfdy(2,1) + CTdminv(1)*(s% e2Phi_bar(2)*dLpdlnT(1) - s% e2Phi_bar(1)*s% L(1)*s% dlnLsdlnT)
        ! core point
        dfdy(2,n) = 0.0
        
        ! J(k+1,k) : 1 <= k <= n-1
        dfdy(3,1:n-1) = -CTdminv(2:n)*s% e2Phi_bar(2:n) * dLpdlnT(1:n-1)
        dfdy(3,n-1) = 0.0   ! for the isothermal core
    end subroutine get_jacobian
    
    subroutine get_num_jacobian(n, x, h, y, f, dfdy, ldfy, lrpar, rpar, lipar, ipar, ierr)
        use const_def, only: dp
        integer, intent(in) :: n, ldfy, lrpar, lipar
        real(dp), intent(in) :: x, h
        real(dp), intent(inout) :: y(n)
        real(dp), intent(out) :: f(n) ! dy/dx
        real(dp), intent(out) :: dfdy(ldfy, n)
        ! dense: dfdy(i, j) = partial f(i) / partial y(j)
        ! banded: dfdy(i-j+mujac+1, j) = partial f(i) / partial y(j)
           ! uses rows 1 to mljac+mujac+1 of dfdy.
           ! The j-th column of the square matrix is stored in the j-th column of the
           ! array dfdy as follows:
           ! dfdy(mujac+1+i-j, j) = partial f(i) / partial y(j)
           ! for max(1, j-mujac)<=i<=min(N, j+mljac)
        integer, intent(inout), pointer :: ipar(:) ! (lipar)
        real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
        integer, intent(out) :: ierr ! nonzero means terminate integration
        type(NScool_info), pointer :: s
        real(dp) :: invtwoh
        ! work arrays
        real(dp), dimension(n) :: ym(n), yp(n), y0(n), fm(n), fp(n)
        integer :: j
        
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
        
        dfdy = 0.0_dp
        invtwoh = 0.5_dp/h
        call get_derivatives(n, x, h, y, f, lrpar, rpar, lipar, ipar, ierr)

        ! first point
        y0 = y
        ym = y0
        yp = y0
        ym(1) = y0(1)-h
        yp(1) = y0(1)+h
        call get_derivatives(n, x, h, ym, fm, lrpar, rpar, lipar, ipar, ierr)
        call get_derivatives(n, x, h, yp, fp, lrpar, rpar, lipar, ipar, ierr)
        dfdy(1,1) = 0.0
        dfdy(2,1) = (fp(1) - fm(1))*invtwoh
        dfdy(3,1) = (fp(2) - fm(2))*invtwoh
        ! interior points
        do j = 2, n-1
           ym = y0
           yp = y0
           ym(j) = y0(j)-h
           yp(j) = y0(j)+h
           call get_derivatives(n, x, h, ym, fm, lrpar, rpar, lipar, ipar, ierr)
           call get_derivatives(n, x, h, yp, fp, lrpar, rpar, lipar, ipar, ierr)
           dfdy(1,j) = (fp(j-1) - fm(j-1))*invtwoh
           dfdy(2,j) = (fp(j) - fm(j))*invtwoh
           dfdy(3,j) = (fp(j+1) - fm(j+1))*invtwoh
        end do
        ! last point
        ym = y0
        yp = y0
        ym(n) = y0(n)-h
        yp(n) = y0(n)+h
        call get_derivatives(n, x, h, ym, fm, lrpar, rpar, lipar, ipar, ierr)
        call get_derivatives(n, x, h, yp, fp, lrpar, rpar, lipar, ipar, ierr)
        dfdy(1,n) = (fp(n-1) - fm(n-1))*invtwoh
        dfdy(2,n) = (fp(n) - fm(n))*invtwoh
        dfdy(3,n) = 0.0
    end subroutine get_num_jacobian
    
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
            
            call interp_value_and_slope(s% tab_lnT, s% n_tab, lnEnu_interp, s% lnT(iz), s% lnenu(iz), s% dlnenu_dlnT(iz), ierr)
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