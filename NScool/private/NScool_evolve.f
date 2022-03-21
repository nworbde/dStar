module NScool_evolve
    use math_lib
    use NScool_def
    
    integer, parameter :: i_id = 1
    integer, parameter :: i_num_terminal_writes = 2
    integer, parameter :: num_deriv_ipar = 2
    integer, parameter :: num_deriv_rpar = 0
    
contains
    subroutine do_integrate_crust(NScool_id,ierr)
        use iso_fortran_env, only : error_unit
        use exceptions_lib
        use constants_def
        use num_def
        use num_lib
        use mtx_def
        use mtx_lib
      
        integer, intent(in) :: NScool_id
        integer, intent(out) :: ierr
        type(NScool_info), pointer :: s
        real(dp) :: t, tend, h, max_step_size
        integer :: which_solver
        integer :: n, max_steps, ijac, itol, nzmax, isparse, mljac, mujac
        integer :: iout, lout, idid
        integer :: lwork, liwork, lrd, lid
        integer :: imas,mumas,mlmas
        real(dp), dimension(1) :: rtol, atol
        real(dp), pointer, dimension(:) :: rpar
        integer, pointer, dimension(:) :: ipar
        real(dp), pointer, dimension(:) :: rpar_decsol, work
        integer, pointer, dimension(:) :: ipar_decsol, iwork
        real(dp), pointer, dimension(:) :: z
        real(dp) :: zmin, zmax
        ! for the bicyclic routine, not used
        real(dp), dimension(:), pointer :: lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk
        integer :: caller_id, nvar, nz
        ! error codes
        character(len=80), dimension(-8:2), parameter :: error_codes = [character(len=80) ::  &
        &   "cannot satisfy given tolerances even after reducing stepsize by 1d30", &
        &   "illegal arg for isolve", &
        &   "terminated by fcn returning nonzero ierr", &
        &   "terminated by jac returning nonzero ierr", &
        &   "matrix is repeatedly singular", &
        &   "step size becomes too small", &
        &   "reached max allowed number of steps", &
        &   "input is not consistent", &
        &   "", &
        &   "computation successful", &
        &   "computation successful (terminated by solout)" ]
        type(alert) :: integration=alert(scope='do_integrate_crust')

        ierr = 0

        call get_NScool_info_ptr(NScool_id,s,ierr)
        if (ierr /= 0) return

        ! set up the equation parameters
        n = s% nz
        t = 0.0_dp
        tend = s% epoch_duration
        h = s% dt
        max_step_size = s% maximum_timestep
        max_steps = s% maximum_number_of_models

        allocate(z(n))
        z = s% lnT(1:s% nz)
        if (s% fix_atmosphere_temperature_when_accreting .and. s% Mdot > 0.0_dp) then
            z(1) = log(s% atmosphere_temperature_when_accreting)
        end if
        
        itol = 0
        ! use coldest temperature to set abolute tolerance
        rtol = s% integration_tolerance
        atol = s% integration_tolerance * minval(z)
        
        ! bounds for temperatures
        zmin = s% min_lg_temperature_integration*ln10
        zmax = s% max_lg_temperature_integration*ln10

        ijac = 1
        nzmax = 0
        isparse = 0
        mljac = 1
        mujac = 1

        imas = 0
        mlmas = 0
        mumas = 0

        caller_id = 0
        nvar = -1
        nz = 0

        iout = 1

        call lapack_work_sizes(n, lrd, lid)
        allocate(rpar_decsol(lrd),ipar_decsol(lid))

        call isolve_work_sizes(n, nzmax, imas, mljac, mujac, mlmas, mumas, liwork, lwork)
        allocate(work(lwork), iwork(liwork))
        iwork = 0
        work = 0.0_dp

        allocate(ipar(num_deriv_ipar),rpar(num_deriv_rpar))

        ipar(i_id) = NScool_id
        ipar(i_num_terminal_writes) = 0
        which_solver = solver_option(trim(s% which_solver),ierr)
        call isolve(which_solver,n,get_derivatives,t,z,tend,h,max_step_size,max_steps, &
          & rtol,atol,itol,zmin,zmax,get_jacobian,ijac,null_sjac,nzmax,isparse, mljac, mujac, &
          & null_mas, imas, mlmas, mumas, evaluate_timestep, iout, lapack_decsol, null_decsols, &
          & null_decsolblk, lrd, rpar_decsol, lid, ipar_decsol,  &
          & caller_id, nvar, nz, lblk, dblk, ublk, uf_lblk, uf_dblk, uf_ublk, null_fcn_blk_dble, null_jac_blk_dble, &
          & work, lwork, iwork, liwork, &
          & num_deriv_rpar, rpar, num_deriv_ipar, ipar, error_unit, idid)

        ! post-mortem
        call integration% report(error_codes(idid))
        if (idid < 0) ierr = idid

        ! print statistics
        write(error_unit,'(a30," = ",i0)') 'num. fcn. evals', iwork(14)
        write(error_unit,'(a30," = ",i0)') 'num. jac. evals', iwork(15)
        write(error_unit,'(a30," = ",i0)') 'num. computed steps', iwork(16)
        write(error_unit,'(a30," = ",i0)') 'num. accepted steps', iwork(17)
        write(error_unit,'(a30," = ",i0)') 'num. rejected steps', iwork(18)
        write(error_unit,'(a30," = ",i0)') 'num. LU decomps.', iwork(19)
        write(error_unit,'(a30," = ",i0)') 'num. forward-backward subs.', iwork(20)

        deallocate(z, rpar_decsol,ipar_decsol,work,iwork)
        nullify(ipar)
        nullify(rpar)

    end subroutine do_integrate_crust

    subroutine evaluate_timestep(nr, xold, x, n, y, rwork_y, iwork_y, interp_y, lrpar, rpar, lipar, ipar, irtrn)
       use iso_fortran_env, only : output_unit, error_unit
       use exceptions_lib
       use NScool_terminal, only : do_write_terminal
       use NScool_history, only : do_write_history
       use NScool_profile, only : do_write_profile
       integer, intent(in) :: nr, n, lrpar, lipar
       real(dp), intent(in) :: xold, x
       real(dp), intent(inout) :: y(:)
       ! y can be modified if necessary to keep it in valid range of possible solutions.
       real(dp), intent(inout), target :: rwork_y(*)
       integer, intent(inout), target :: iwork_y(*)
       integer, intent(inout), pointer :: ipar(:) ! (lipar)
       real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
       interface
          include 'num_interp_y.dek'
       end interface
       integer, intent(out) :: irtrn ! < 0 causes solver to return to calling program.
       type(NScool_info), pointer :: s
       integer :: ierr
       logical :: print_terminal_header
       character(len=256) :: filename      
       type(assertion) :: got_pointer=assertion( &
       &   scope='evaluate_timestep',message='unable to access NScool_info')
       type(assertion) :: equations=assertion( &
       &   scope='evaluate_timestep',message='wrong number of equations')
       type(failure) :: coefficient_error=failure(scope='evaluate_timestep', &
       &    message='evaluatning coefficients')
       type(failure) :: luminosity_error=failure(scope='evaluate_timestep', &
       &    message='evaluating luminosity')
       type(warning) :: terminal_warning=warning(scope='evaluate_timestep', &
       &    message='writing to terminal')
       type(warning) :: history_warning=warning(scope='evaluate_timestep', &
       &    message='writing history')
       type(warning) :: profile_warning=warning(scope='evaluate_timestep', &
       &    message='writing profile')
       type(alert) :: writing=alert(scope='evaluate_timestep')
       character(len=128) :: status_message
           
       irtrn = 0
       ierr = 0

       ! get the crust information information
       call get_NScool_info_ptr(ipar(i_id), s, ierr)
       call got_pointer% assert(ierr == 0)
       call equations% assert(n == s% nz)
       
       if (s% suppress_first_step_output .and. x == xold) return

       s% model = nr + s% starting_number_for_profile
       s% tsec = x + s% epoch_start_time
       s% dt = x-xold
  
       s% lnT(1:n) = y(1:n)
       s% T(1:n) = exp(s% lnT(1:n))
       call interpolate_temps(s)
  
       call get_coefficients(s,ierr)
       if (coefficient_error% raised(ierr)) then
           irtrn=-1
           return
       end if
  
       call evaluate_luminosity(s, ierr)
       if (luminosity_error% raised(ierr)) then
           irtrn=-1
           return
       end if

       ! update terminal information
       if (mod(s% model,s% write_interval_for_terminal) == 0) then
          if (mod(ipar(i_num_terminal_writes),s% write_interval_for_terminal_header) == 0) then
             print_terminal_header = .TRUE.
          else
             print_terminal_header = .FALSE.
          end if
          call do_write_terminal(ipar(i_id), ierr, print_terminal_header)
          if (terminal_warning% raised(ierr)) return
          ipar(i_num_terminal_writes) = ipar(i_num_terminal_writes) + 1
       end if
  
       ! update history log
       if (mod(s% model, s% write_interval_for_history) == 0) then
          call do_write_history(ipar(i_id),ierr)
          if (history_warning% raised(ierr)) return
          write (status_message,'(a,i0,a)') 'saving model ',s% model,' to history log'
          call writing% report(status_message)
       end if
  
       ! update profile log
       if (mod(s% model, s% write_interval_for_profile) == 0) then
          write(filename,'(a,i4.4)') 'profile',s% model
          call do_write_profile(ipar(i_id),ierr)
          if (profile_warning% raised(ierr)) return
          write (status_message,'(a,i0,a)') 'saving model ',s% model,' to profile log'
          call writing% report(status_message)
       end if
    end subroutine evaluate_timestep

    subroutine get_derivatives(n, x, h, y, f, lrpar, rpar, lipar, ipar, ierr)
        use const_def, only: dp
        use exceptions_lib
        integer, intent(in) :: n, lrpar, lipar
        real(dp), intent(in) :: x, h
        real(dp), intent(inout) :: y(:) ! okay to edit y if necessary (e.g., replace negative values by zeros)
        real(dp), intent(inout) :: f(:) ! dy/dx
        integer, intent(inout), pointer :: ipar(:) ! (lipar)
        real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
        integer, intent(out) :: ierr ! nonzero means retry with smaller timestep.
        type(NScool_info), pointer :: s
        type(assertion) :: got_pointer=assertion( &
        &   scope='get_derivatives',message='unable to access NScool_info')
        type(assertion) :: equations=assertion( &
        &   scope='get_derivatives',message='wrong number of equations')

        ierr = 0
        call get_NScool_info_ptr(ipar(i_id), s, ierr)
        ! fatal erors
        call got_pointer% assert(ierr == 0)
        call equations% assert(n == s% nz)

        s% lnT(1:n) = y(1:n)
        s% T(1:n) = exp(s% lnT(1:n))

        call interpolate_temps(s)

        call get_coefficients(s, ierr)
        call evaluate_luminosity(s,ierr)

        f(1:n-1) = ((s% L(2:n)*s% e2Phi_bar(2:n) - s% L(1:n-1)*s% e2Phi_bar(1:n-1))/s% dm(1:n-1)/s% e2Phi(1:n-1)  &
        &   + s% enuc(1:n-1) - s% enu(1:n-1))*s% ePhi(1:n-1)/s% Cp(1:n-1)/s% T(1:n-1)
        if (s% make_inner_boundary_insulating .and. .not. s% fix_core_temperature) then
            f(n) = (( -s% L(n)*s% e2Phi_bar(n))/s% dm(n)/s% e2Phi(n)  &
            &   + s% enuc(n) - s% enu(n))*s% ePhi(n)/s% Cp(n)/s% T(n)
        else
            f(n) = 0.0_dp
         end if
        
        if (s% fix_atmosphere_temperature_when_accreting .and. s% Mdot > 0.0_dp) then
            f(1) = 0.0_dp
        end if
        
    end subroutine get_derivatives
    
    subroutine get_jacobian(n, x, h, y, f, dfdy, ldfy, lrpar, rpar, lipar, ipar, ierr)
        use const_def, only: dp
        use exceptions_lib
        integer, intent(in) :: n, ldfy, lrpar, lipar
        real(dp), intent(in) :: x, h
        real(dp), intent(inout) :: y(:)
        real(dp), intent(inout) :: f(:) ! dy/dx
        real(dp), intent(inout) :: dfdy(:,:) !dfdy(ldfy, n)
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
        type(assertion) :: got_pointer=assertion( &
        &   scope='get_jacobian',message='unable to access NScool_info')
        type(assertion) :: equations=assertion( &
        &   scope='get_jacobian',message='wrong number of equations')
        
        ierr = 0
        call get_NScool_info_ptr(ipar(i_id), s, ierr)
        ! fatal erors
        call got_pointer% assert(ierr == 0)
        call equations% assert(n == s% nz)
        
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
        if (s% make_inner_boundary_insulating .and. .not. s% fix_core_temperature) then
            dfdy(2,n) = dfdy(2,n) + CTdminv(n)*(- s% e2Phi_bar(n)*dLdlnT(n))
        else
            dfdy(2,n) = 0.0
        end if
        
        ! J(k+1,k) : 1 <= k <= n-1
        dfdy(3,1:n-1) = -CTdminv(2:n)*s% e2Phi_bar(2:n) * dLpdlnT(1:n-1)
        dfdy(3,n-1) = 0.0   ! for the isothermal core

        if (s% fix_atmosphere_temperature_when_accreting .and. s% Mdot > 0.0_dp) then
            dfdy(1,2) = 0.0_dp
            dfdy(2:3,1) = 0.0_dp
        end if
    end subroutine get_jacobian
    
    subroutine get_num_jacobian(n, x, h, y, f, dfdy, ldfy, lrpar, rpar, lipar, ipar, ierr)
        use const_def, only: dp
        use exceptions_lib
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
        type(assertion) :: got_pointer=assertion( &
        &   scope='get_num_jacobian',message='unable to access NScool_info')
        type(assertion) :: equations=assertion( &
        &   scope='get_num_jacobian',message='wrong number of equations')
        
        ierr = 0
        call get_NScool_info_ptr(ipar(i_id), s, ierr)
        call got_pointer% assert(ierr == 0)
        call equations% assert(n == s% nz)
        
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
        s% Lsurf = s% area(1) * exp10(lgflux)     ! emergent luminosity
        s% dlnLsdlnT = dlgflux
        s% Teff = exp10(lgTeff)     ! surface effective temperature
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
        use exceptions_lib
        type(NScool_info), pointer :: s
        integer, intent(out) :: ierr
        integer :: iz
        real(dp), dimension(:), pointer :: lnKcond_interp, lnCp_interp, lnGamma_interp, lnEnu_interp
        type(failure) :: interpolation_error=failure(scope='get_coefficients')
        type(failure) :: heating_error=failure(scope='get_coefficients', &
        &   message='calculating nuclear heating')
        
        do iz = 1, s% nz
            
            lnKcond_interp(1:4*s% n_tab) => s% tab_lnK(1:4*s% n_tab, iz)
            lnCp_interp(1:4*s% n_tab) => s% tab_lnCp(1:4*s% n_tab, iz)
            lnGamma_interp(1:4*s% n_tab) => s% tab_lnGamma(1:4*s% n_tab, iz)
            lnEnu_interp(1:4*s% n_tab) => s% tab_lnEnu(1:4*s% n_tab, iz)
        
            call interp_value_and_slope(s% tab_lnT, s% n_tab, lnKcond_interp, s% lnT_bar(iz), s% lnK(iz), s% dlnK_dlnT(iz), ierr)
            if (interpolation_error% raised(ierr)) then
                call report('lnK', iz)
                return
            end if

            call interp_value_and_slope(s% tab_lnT, s% n_tab, lnCp_interp, s% lnT(iz), s% lnCp(iz), s% dlnCp_dlnT(iz), ierr)
            if (interpolation_error% raised(ierr)) then
                call report('lnCp',iz)
                return
            end if
            
            call interp_value_and_slope(s% tab_lnT, s% n_tab, lnGamma_interp, s% lnT(iz), s% lnGamma(iz),  &
            &   s% dlnGamma_dlnT(iz), ierr)
            if (interpolation_error% raised(ierr)) then
                call report('lnGamma',iz)
                return
            end if
            
            call interp_value_and_slope(s% tab_lnT, s% n_tab, lnEnu_interp, s% lnT(iz), s% lnenu(iz), s% dlnenu_dlnT(iz), ierr)
            if (interpolation_error% raised(ierr)) then
                call report('lnEnu',iz)
                return
            end if
            
            s% Kcond = exp(s% lnK)
            s% Cp = exp(s% lnCp)
            s% Gamma = exp(s% lnGamma)
            s% enu = exp(s% lnenu)
            
        end do
        
        call get_nuclear_heating(s, ierr)
        if (heating_error% raised(ierr)) return
        
        contains
        subroutine report(str,izone)
           character(len=*), intent(in) :: str
           integer, intent(in) :: izone
           type(alert) :: status=alert(scope='get_coefficient',level=0)
           character(len=128) :: msg
           write (msg,'(a,i0)') 'interpolating '//str//', zone ',iz
           call status% report(msg)
        end subroutine report
            
    end subroutine get_coefficients
    
    subroutine get_nuclear_heating(s, ierr)
        use constants_def
        type(NScool_info), pointer :: s
        integer, intent(out) :: ierr
        real(dp) :: Ptop, Pbot, Q
        
        ierr = 0
        s% enuc = 0.0_dp
        
        ! outer crust
        Ptop = exp10(s% lgP_min_heating_outer)
        Pbot = exp10(s% lgP_max_heating_outer)
        Q = s% Q_heating_outer
        call do_one
        
        ! inner crust
        Ptop = exp10(s% lgP_min_heating_inner)
        Pbot = exp10(s% lgP_max_heating_inner)
        Q = s% Q_heating_inner
        call do_one
        
        ! extra, if it exists
        if (s% turn_on_extra_heating) then
            Ptop = exp10(s% lgP_min_heating_shallow)
            Pbot = exp10(s% lgP_max_heating_shallow)
            Q = s% Q_heating_shallow
            call do_one
        end if
        
        ! now hook in user-supplied routine
        if (s% use_other_set_heating) then
            call s% other_set_heating(s% id, ierr)
        end if
        
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
