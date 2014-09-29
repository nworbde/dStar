program test_NScool
    use constants_def
    use NScool_def
    use NScool_lib
    use create_model
    use NScool_ctrls_io, only: write_controls
    use NScool_evolve
    use NScool_profile
    use, intrinsic :: iso_fortran_env, only: output_unit
    use interp_1d_lib
    
    character(len=*), parameter :: my_dStar_dir = '../../../dStar'
    character(len=*), parameter :: inlist = 'test_inlist'
    type(NScool_info), pointer :: s
    integer :: ierr, NScool_id,i
    real(dp), dimension(:,:), pointer :: lnEnu_val, lnKcond_val, lnCp_val
    real(dp), dimension(:), pointer :: lnEnu_interp, lnKcond_interp, lnCp_interp
    real(dp) :: lnT, lnK, lnEnu, lnC, dlnK, dlnEnu, dlnC
    real(dp), dimension(:), pointer :: rpar
    integer, dimension(:), pointer :: ipar
    real(dp) :: h
    real(dp), dimension(:), allocatable :: y, f
    
    call NScool_init(my_dStar_dir, ierr)
    call check_okay('NScool_init',ierr)
    
    NScool_id = alloc_NScool(ierr)
    call check_okay('NScool_id',ierr)
    
    call NScool_setup(NScool_id,inlist,ierr)
    call check_okay('NScool_setup',ierr)
    ierr = 0
!     call write_controls(output_unit,ierr)
    
    call get_NScool_info_ptr(NScool_id,s,ierr)
    call do_setup_crust_zones(s, ierr)
    call check_okay('do_setup_crust_zones',ierr)
    call do_setup_crust_composition(s, ierr)
    call check_okay('do_setup_crust_composition',ierr)
    call do_setup_crust_transport(s, ierr)
    call check_okay('do_setup_crust_transport',ierr)
    
    s% Mdot = 1.0e16
    call get_coefficients(s, ierr)
    call check_okay('get_coefficients',ierr)
    call interpolate_temps(s)
    call evaluate_luminosity(s,ierr)
    call check_okay('evaluate_luminosity',ierr)
    
    call do_write_profile(NScool_id,ierr)
    call check_okay('do_write_profile',ierr)
    
    h = 0.1
    allocate(y(s% nz), f(s% nz))
    y = log(s% T)
    allocate(rpar(num_deriv_rpar), ipar(num_deriv_ipar))
    ipar(i_id) = NScool_id
    call get_derivatives(s% nz, 0.0_dp, h, y, f, num_deriv_rpar, rpar, num_deriv_ipar, ipar, ierr)

    do i = 1, s% nz
        write (output_unit,'(i5,es16.8)') i, 1.0_dp/f(i)
!         lnCp_val(1:4,1:s% n_tab) => s% tab_lnCp(1:4*s% n_tab, i)
!         lnEnu_val(1:4,1:s% n_tab) => s% tab_lnEnu(1:4*s% n_tab, i)
!         lnKcond_val(1:4,1:s% n_tab) => s% tab_lnK(1:4*s% n_tab, i)
!         lnKcond_interp(1:4*s% n_tab) => s% tab_lnK(1:4*s% n_tab, i)
!         lnCp_interp(1:4*s% n_tab) => s% tab_lnCp(1:4*s% n_tab, i)
!         lnEnu_interp(1:4*s% n_tab) => s% tab_lnEnu(1:4*s% n_tab, i)
!
!         lnT = log(s% T_bar(i))
!         call interp_value_and_slope(s% tab_lnT, s% n_tab, lnKcond_interp, lnT, lnK, dlnK, ierr)
!         if (ierr /= 0) then
!             print *, 'bad interp in lnK'
!         end if
!         lnT = log(s% T(i))
!         call interp_value_and_slope(s% tab_lnT, s% n_tab, lnCp_interp, lnT, lnC, dlnC, ierr)
!         if (ierr /= 0) then
!             print *, 'bad interp in lnC'
!         end if
!         call interp_value_and_slope(s% tab_lnT, s% n_tab, lnEnu_interp, lnT, lnEnu, dlnEnu, ierr)
!         if (ierr /= 0) then
!             print *, 'bad interp in lnEnu'
!         end if
!         write (output_unit,'(4es15.8,3(f12.8))') s% dm_bar(i), s% rho_bar(i), s% P_bar(i), s% area(i), &
!         &    lnKcond_val(1,42)/ln10, lnK/ln10, dlnK
!         write (output_unit,'(t8,3es15.8,2(f14.10),2(f5.1),7(f14.10),es15.8)') s% dm(i), s% rho(i), s% T(i), s% ePhi(i), &
!         &   s% eLambda(i), &
!         &   s% ionic(i)% Z, s% ionic(i)% A, s% tab_lnT(42)/ln10, lnCp_val(1,42)/ln10, lnC/ln10, dlnC, lnEnu_val(1,42)/ln10, &
!         &   lnEnu/ln10, dlnEnu, s% enuc(i)
! !         write (output_unit,'(2es15.8,tr2,19(f10.6))') s% P_bar(i), s% T_bar(i), s% Yion_bar(1:s% ncharged,i),s% Xneut_bar(i)
! !         write (output_unit,'(2es15.8,tr2,f10.6,)') s% P_bar(i), s% T_bar(i), s% Xneut_bar(i)
    end do
! !     write (output_unit,'(2es15.8,tr2,19(f10.6))') s% P_bar(i), s% T_bar(i), s% Yion_bar(1:s% ncharged,i), s% Xneut_bar(i)
! ! !     write (output_unit,'(es15.8,tr2,f14.10,tr1,es15.8,f14.10)') s% P_bar(s% nz), s% ePhi_bar(s% nz), s% m(s% nz),  &
! ! !     &   s% eLambda_bar(s% nz)
! !
! !     write (output_unit,*) s% charged_ids
    
!
!     write (output_unit, '(/,/,"L = ",es15.8,f12.8)') dot_product(s% dm, s% enuc),  &
!     &   dot_product(s% dm, s% enuc) * ergs_to_mev/ s% Mdot /avogadro
    
    deallocate(rpar, ipar, y, f)
    
    call NScool_shutdown
    
contains
	subroutine check_okay(msg,ierr)
		use iso_fortran_env, only : error_unit
		character(len=*), intent(in) :: msg
		integer, intent(inout) :: ierr
		if (ierr /= 0) then
			write (error_unit,*) trim(msg)//': ierr = ',ierr
			if (ierr < 0) stop
		end if
	end subroutine check_okay
end program test_NScool
