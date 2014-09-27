program test_NScool
    use NScool_def
    use NScool_lib
    use create_model
    use NScool_ctrls_io, only: write_controls
    use, intrinsic :: iso_fortran_env, only: output_unit
    
    character(len=*), parameter :: my_dStar_dir = '../../../dStar'
    character(len=*), parameter :: inlist = 'test_inlist'
    type(NScool_info), pointer :: s
    integer :: ierr, NScool_id,i
    real(dp), dimension(:,:), pointer :: lnEnu_val, lnKcond_val, lnCp_val
    
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
    
!     call tov_write_crust

    do i = 1, s% nz
        lnCp_val(1:4,1:s% n_tab) => s% tab_lnCp(1:4*s% n_tab, i)
        lnEnu_val(1:4,1:s% n_tab) => s% tab_lnEnu(1:4*s% n_tab, i)
        lnKcond_val(1:4,1:s% n_tab) => s% tab_lnK(1:4*s% n_tab, i)

        write (output_unit,'(3es15.8,f12.8)') s% dm_bar(i), s% rho_bar(i), s% P_bar(i), lnKcond_val(1,42)/ln10
        write (output_unit,'(t8,3es15.8,2(f14.10),2(f5.1),3(f14.10))') s% dm(i), s% rho(i), s% T(i), s% ePhi(i), &
        &   s% eLambda(i), &
        &   s% ionic(i)% Z, s% ionic(i)% A, s% tab_lnT(42)/ln10, lnCp_val(1,42)/ln10, lnEnu_val(1,42)/ln10
!         write (output_unit,'(2es15.8,tr2,19(f10.6))') s% P_bar(i), s% T_bar(i), s% Yion_bar(1:s% ncharged,i),s% Xneut_bar(i)
!         write (output_unit,'(2es15.8,tr2,f10.6,)') s% P_bar(i), s% T_bar(i), s% Xneut_bar(i)
    end do
!     write (output_unit,'(2es15.8,tr2,19(f10.6))') s% P_bar(i), s% T_bar(i), s% Yion_bar(1:s% ncharged,i), s% Xneut_bar(i)
! !     write (output_unit,'(es15.8,tr2,f14.10,tr1,es15.8,f14.10)') s% P_bar(s% nz), s% ePhi_bar(s% nz), s% m(s% nz),  &
! !     &   s% eLambda_bar(s% nz)
!
!     write (output_unit,*) s% charged_ids
    
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
