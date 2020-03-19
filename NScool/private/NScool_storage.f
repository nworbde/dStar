module storage
    use NScool_private_def

    integer, parameter :: do_deallocate = 0
    integer, parameter :: do_allocate = 1
    
contains
    subroutine free_NScool_data(id, ierr)
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        type(NScool_info), pointer :: s
        call get_NScool_info_ptr(id,s,ierr)
        if (ierr /= 0) return
        call free_arrays(s)
    end subroutine free_NScool_data

    subroutine free_arrays(s)
        type(NScool_info), pointer :: s
        integer :: ierr
        
        call NScool_iso_arrays(s, do_deallocate, ierr)
        call NScool_info_arrays(s, do_deallocate, ierr)
        call NScool_work_arrays(s, do_deallocate, ierr)
        call NScool_epoch_arrays(s, do_deallocate, ierr)
        call free_NScool(s)
    end subroutine free_arrays
    
    subroutine allocate_NScool_iso_arrays(s, ierr)
        type(NScool_info), pointer :: s
        integer, intent(out) :: ierr
        call NScool_iso_arrays(s, do_allocate, ierr)
    end subroutine allocate_NScool_iso_arrays

    subroutine free_NScool_iso_arrays(s, ierr)
        type(NScool_info), pointer :: s
        integer, intent(out) :: ierr
        call NScool_info_arrays(s, do_deallocate, ierr)
    end subroutine free_NScool_iso_arrays

    subroutine allocate_NScool_info_arrays(s, ierr)
        type(NScool_info), pointer :: s
        integer, intent(out) :: ierr
        call NScool_info_arrays(s, do_allocate, ierr)
    end subroutine allocate_NScool_info_arrays

    subroutine free_NScool_info_arrays(s, ierr)
        type(NScool_info), pointer :: s
        integer, intent(out) :: ierr
        call NScool_info_arrays(s, do_deallocate, ierr)
    end subroutine free_NScool_info_arrays

    subroutine allocate_NScool_work_arrays(s, ierr)
        type(NScool_info), pointer :: s
        integer, intent(out) :: ierr
        call NScool_work_arrays(s, do_allocate, ierr)
    end subroutine allocate_NScool_work_arrays

    subroutine free_NScool_work_arrays(s, ierr)
        type(NScool_info), pointer :: s
        integer, intent(out) :: ierr
        call NScool_work_arrays(s, do_deallocate, ierr)
    end subroutine free_NScool_work_arrays
    
    subroutine allocate_NScool_epoch_arrays(s, ierr)
        type(NScool_info), pointer :: s
        integer, intent(out) :: ierr
        call NScool_epoch_arrays(s, do_allocate, ierr)
    end subroutine allocate_NScool_epoch_arrays

    subroutine free_NScool_epoch_arrays(s, ierr)
        type(NScool_info), pointer :: s
        integer, intent(out) :: ierr
        call NScool_epoch_arrays(s, do_deallocate, ierr)
    end subroutine free_NScool_epoch_arrays
        
    subroutine NScool_iso_arrays(s, action, ierr)
        use exceptions_lib
        type(NScool_info), pointer :: s
        integer, intent(in) :: action
        integer, intent(out) :: ierr
        integer :: nisos
        type(failure) :: iso_failure=failure(scope='NScool_iso_arrays')
        ierr = 0
        nisos = s% nisos
        do
!             call do1(s% iso_ids)
!             if (failed('iso_ids')) exit
            call do1(s% charged_ids)
            if (iso_failure% raised(ierr,'charged_ids')) exit
            return
        end do
        ierr = -1

    contains
        subroutine do1(ptr)
            integer, dimension(:), pointer :: ptr
            call do1Dint(ptr, nisos, action, ierr)
        end subroutine do1
    end subroutine NScool_iso_arrays

    subroutine NScool_info_arrays(s, action, ierr)
        use exceptions_lib
        type(NScool_info), pointer :: s
        integer, intent(in) :: action
        integer, intent(out) :: ierr
        integer :: nz, nisos !, nepochs
        type(failure) :: info_failure=failure(scope='NScool_info_arrays')
        
        ierr = 0
        nz = s% nz
        nisos = s% nisos
!         nepochs = s% number_epochs
        
        do
            call do1(s% dm)
            if (info_failure% raised(ierr,'dm')) exit
!             call do2(s% X)
!             if (info_failure% raised(ierr,'X')) exit
            call do2(s% Yion)
            if (info_failure% raised(ierr,'Yion')) exit
            call do1(s% Xneut)
            if (info_failure% raised(ierr,'Xneut')) exit
            call do1c(s% ionic)
            if (info_failure% raised(ierr,'ionic')) exit
            call do1(s% P)
            if (info_failure% raised(ierr,'P')) exit
            call do1(s% lnT)
            if (info_failure% raised(ierr,'lnT')) exit
            call do1(s% T)
            if (info_failure% raised(ierr,'T')) exit
            call do1(s% ePhi)
            if (info_failure% raised(ierr,'ePhi')) exit
            call do1(s% e2Phi)
            if (info_failure% raised(ierr,'e2Phi')) exit
            call do1(s% eLambda)
            if (info_failure% raised(ierr,'eLambda')) exit
            call do1(s% rho)
            if (info_failure% raised(ierr,'rho')) exit
            call do1(s% lnCp)
            if (info_failure% raised(ierr,'lnCp')) exit
            call do1(s% Cp)
            if (info_failure% raised(ierr,'Cp')) exit
            call do1(s% dlnCp_dlnT)
            if (info_failure% raised(ierr,'dlnCp_dlnT')) exit
            call do1(s% lnGamma)
            if (info_failure% raised(ierr,'lnGamma')) exit
            call do1(s% Gamma)
            if (info_failure% raised(ierr,'Gamma')) exit
            call do1(s% dlnGamma_dlnT)
            if (info_failure% raised(ierr,'dlnGamma_dlnT')) exit
            call do1(s% enu)
            if (info_failure% raised(ierr,'enu')) exit
            call do1(s% lnenu)
            if (info_failure% raised(ierr,'lnenu')) exit
            call do1(s% dlnenu_dlnT)
            if (info_failure% raised(ierr,'dlnenu_dlnT')) exit
            call do1(s% enuc)
            if (info_failure% raised(ierr,'enuc')) exit
            
            call do1(s% m)
            if (info_failure% raised(ierr,'m')) exit
            call do1(s% area)
            if (info_failure% raised(ierr,'area')) exit
            call do1(s% L)
            if (info_failure% raised(ierr,'L')) exit
            call do1(s% dm_bar)
            if (info_failure% raised(ierr,'dm_bar')) exit
!             call do2(s% X_bar)
!             if (failed('X_bar')) exit
            call do2(s% Yion_bar)
            if (info_failure% raised(ierr,'Yion_bar')) exit
            call do1(s% Xneut_bar)
            if (info_failure% raised(ierr,'Xneut_bar')) exit
            call do1c(s% ionic_bar)
            if (info_failure% raised(ierr,'ionic_bar')) exit
            call do1(s% P_bar)
            if (info_failure% raised(ierr,'P_bar')) exit
            call do1(s% lnT_bar)
            if (info_failure% raised(ierr,'lnT_bar')) exit
            call do1(s% T_bar)
            if (info_failure% raised(ierr,'T_bar')) exit
            call do1(s% ePhi_bar)
            if (info_failure% raised(ierr,'ePhi_bar')) exit
            call do1(s% e2Phi_bar)
            if (info_failure% raised(ierr,'e2Phi_bar')) exit
            call do1(s% eLambda_bar)
            if (info_failure% raised(ierr,'eLambda_bar')) exit
            call do1(s% rho_bar)
            if (info_failure% raised(ierr,'rho_bar')) exit
            call do1(s% Kcond)
            if (info_failure% raised(ierr,'Kcond')) exit
            call do1(s% lnK)
            if (info_failure% raised(ierr,'lnK')) exit
            call do1(s% dlnK_dlnT)
            if (info_failure% raised(ierr,'dlnK_dlnT')) exit
!             call do1epoch(s% t_monitor)
!             if (failed('t_monitor')) exit
!             call do1epoch(s% Teff_monitor)
!             if (failed('Teff_monitor')) exit
!             call do1epoch(s% Qb_monitor)
!             if (failed('Qb_monitor')) exit
            return
        end do
        ierr = -1
        
    contains
!         subroutine do1epoch(ptr)
!             real(dp), dimension(:), pointer :: ptr
!             call do1D(ptr, nepochs, action, ierr)
!         end subroutine do1epoch
        subroutine do1(ptr)
            real(dp), dimension(:), pointer :: ptr
            call do1D(ptr, nz, action, ierr)
        end subroutine do1
        subroutine do2(ptr)
            real(dp), dimension(:,:), pointer :: ptr
            call do2D(ptr, nisos, nz, action, ierr)
        end subroutine do2
        subroutine do1c(ptr)
            use nucchem_def, only: composition_info_type
            type(composition_info_type), dimension(:), pointer :: ptr
            call do1Dcomp(ptr,nz,action, ierr)
        end subroutine do1c
    end subroutine NScool_info_arrays

    subroutine NScool_work_arrays(s,action,ierr)
        use exceptions_lib
        type(NScool_info), pointer :: s
        integer, intent(in) :: action
        integer, intent(out) :: ierr
        integer :: ntab, nz
        type(failure) :: work_failure=failure(scope='NScool_work_arrays')
        
        ierr = 0
        nz = s% nz
        ntab = s% n_tab
        do
            call do1(s% tab_lnT)
            if (work_failure% raised(ierr,'tab_lnT')) exit
            call do2(s% tab_lnenu)
            if (work_failure% raised(ierr,'tab_lnenu')) exit
            call do2(s% tab_lnCp)
            if (work_failure% raised(ierr,'tab_lnCp')) exit
            call do2(s% tab_lnGamma)
            if (work_failure% raised(ierr,'tab_lnGamma')) exit
            call do2(s% tab_lnK)
            if (work_failure% raised(ierr,'tab_lnK')) exit
            return
        end do
        ierr = -1
    contains
        subroutine do1(ptr)
            real(dp), dimension(:), pointer :: ptr
            call do1D(ptr,ntab,action,ierr)
        end subroutine do1
        subroutine do2(ptr)
            real(dp), dimension(:,:), pointer :: ptr
            call do2D(ptr,4*ntab,nz,action,ierr)
        end subroutine do2
    end subroutine NScool_work_arrays
    
    subroutine NScool_epoch_arrays(s, action, ierr)
        use exceptions_lib
        type(NScool_info), pointer :: s
        integer, intent(in) :: action
        integer, intent(out) :: ierr
        integer :: nepochs
        type(failure) :: epoch_failure=failure(scope='NScool_epoch_arrays')
        
        ierr = 0
        nepochs = s% number_epochs        
        do
            call do1(s% epoch_Mdots)
            if (epoch_failure% raised(ierr,'epoch_Mdots')) exit
            call do1(s% epoch_boundaries,start=0)
            if (epoch_failure% raised(ierr,'epoch_boundaries')) exit
            call do1(s% t_monitor)
            if (epoch_failure% raised(ierr,'t_monitor')) exit
            call do1(s% Teff_monitor)
            if (epoch_failure% raised(ierr,'Teff_monitor')) exit
            call do1(s% Qb_monitor)
            if (epoch_failure% raised(ierr,'Qb_monitor')) exit
            return
        end do
        ierr = -1
    contains
        subroutine do1(ptr,start)
            real(dp), dimension(:), pointer :: ptr
            integer, optional, intent(in) :: start
            if (present(start)) then
                call do1Doffset(ptr,start,nepochs,action,ierr)
            else
                call do1D(ptr,nepochs,action,ierr)
            end if
        end subroutine do1
    end subroutine NScool_epoch_arrays
    
    subroutine do1D(ptr,sz1,action,ierr)
        real(dp), dimension(:), pointer :: ptr
        integer, intent(in) :: sz1,action
        integer, intent(out) :: ierr
        ierr = 0
        select case(action)
            case (do_deallocate)
            if (associated(ptr)) then
                deallocate(ptr)
                nullify(ptr)
            end if
            case (do_allocate)
            allocate(ptr(sz1),stat=ierr)
        end select
    end subroutine do1D
    
    subroutine do1Doffset(ptr,st1,sz1,action,ierr)
        real(dp), dimension(:), pointer :: ptr
        integer, intent(in) :: st1,sz1,action
        integer, intent(out) :: ierr
        ierr = 0
        select case(action)
            case (do_deallocate)
            if (associated(ptr)) then
                deallocate(ptr)
                nullify(ptr)
            end if
            case (do_allocate)
            allocate(ptr(st1:sz1),stat=ierr)
        end select
    end subroutine do1Doffset
    
    subroutine do1Dint(ptr,sz1,action,ierr)
        integer, dimension(:), pointer :: ptr
        integer, intent(in) :: sz1,action
        integer, intent(out) :: ierr
        ierr = 0
        select case(action)
            case (do_deallocate)
            if (associated(ptr)) then
                deallocate(ptr)
                nullify(ptr)
            end if
            case (do_allocate)
            allocate(ptr(sz1),stat=ierr)
        end select
    end subroutine do1Dint

    subroutine do1Dcomp(ptr,sz1,action,ierr)
        type(composition_info_type), dimension(:), pointer :: ptr
        integer, intent(in) :: sz1,action
        integer, intent(out) :: ierr
        ierr = 0
        select case(action)
            case (do_deallocate)
            if (associated(ptr)) then
                deallocate(ptr)
                nullify(ptr)
            end if
            case (do_allocate)
            allocate(ptr(sz1),stat=ierr)
        end select
    end subroutine do1Dcomp
    
    subroutine do1Depoch(ptr,lb,ub,action,ierr)
        real(dp), dimension(:), pointer :: ptr
        integer, intent(in) :: lb, ub, action
        integer, intent(out) :: ierr
        ierr = 0
        select case(action)
            case (do_deallocate)
            if (associated(ptr)) then
                deallocate(ptr)
                nullify(ptr)
            end if
            case (do_allocate)
            allocate(ptr(lb:ub),stat=ierr)
        end select
    end subroutine do1Depoch
    
    subroutine do2D(ptr,sz1,sz2,action,ierr)
        real(dp), dimension(:,:), pointer :: ptr
        integer, intent(in) :: sz1,sz2,action
        integer, intent(out) :: ierr
        ierr = 0
        select case(action)
            case (do_deallocate)
            if (associated(ptr)) then
                deallocate(ptr)
                nullify(ptr)
            end if
            case (do_allocate)
            allocate(ptr(sz1,sz2),stat=ierr)
        end select
    end subroutine do2D

end module storage
