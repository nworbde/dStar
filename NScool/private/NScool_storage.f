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
        
    subroutine NScool_iso_arrays(s, action, ierr)
        type(NScool_info), pointer :: s
        integer, intent(in) :: action
        integer, intent(out) :: ierr
        integer :: nisos
        
        ierr = 0
        nisos = s% nisos
        do
            call do1(s% iso_ids)
            if (failed('iso_ids')) exit
            call do1(s% charged_ids)
            if (failed('charged_ids')) exit
            return
        end do
        ierr = -1

    contains
        subroutine do1(ptr)
            integer, dimension(:), pointer :: ptr
            call do1Dint(ptr, nisos, action, ierr)
        end subroutine do1
        function failed(str)
            character(len=*), intent(in) :: str
            logical :: failed
            failed = .FALSE.
            if (ierr == 0) return
            write (*,*) 'NScool_info_arrays failed for ',trim(str)
            failed = .TRUE.
        end function failed
    end subroutine NScool_iso_arrays

    subroutine NScool_info_arrays(s, action, ierr)
        type(NScool_info), pointer :: s
        integer, intent(in) :: action
        integer, intent(out) :: ierr
        integer :: nz, nisos
        
        ierr = 0
        nz = s% nz
        nisos = s% nisos
        
        do
            call do1(s% dm)
            if (failed('dm')) exit
            call do2(s% X)
            if (failed('X')) exit
            call do2(s% Yion)
            if (failed('Yion')) exit
            call do1c(s% ionic)
            if (failed('ionic')) exit
            call do1(s% P)
            if (failed('P')) exit
            call do1(s% lnT)
            if (failed('lnT')) exit
            call do1(s% T)
            if (failed('T')) exit
            call do1(s% phi)
            if (failed('phi')) exit
            call do1(s% eLambda)
            if (failed('eLambda')) exit
            call do1(s% rho)
            if (failed('rho')) exit
            call do1(s% lnCp)
            if (failed('lnCp')) exit
            call do1(s% Cp)
            if (failed('Cp')) exit
            call do1(s% dlnCp_dlnT)
            if (failed('dlnCp_dlnT')) exit
            call do1(s% enu)
            if (failed('enu')) exit
            call do1(s% lnenu)
            if (failed('lnenu')) exit
            call do1(s% dlnenu_dlnT)
            if (failed('dlnenu_dlnT')) exit
            call do1(s% enuc)
            if (failed('enuc')) exit
            
            call do1(s% m)
            if (failed('m')) exit
            call do1(s% L)
            if (failed('L')) exit
            call do1(s% dm_bar)
            if (failed('dm_bar')) exit
            call do2(s% X_bar)
            if (failed('X_bar')) exit
            call do2(s% Yion_bar)
            if (failed('Yion_bar')) exit
            call do1c(s% ionic_bar)
            if (failed('ionic_bar')) exit
            call do1(s% P_bar)
            if (failed('P_bar')) exit
            call do1(s% lnT_bar)
            if (failed('lnT_bar')) exit
            call do1(s% T_bar)
            if (failed('T_bar')) exit
            call do1(s% phi_bar)
            if (failed('phi_bar')) exit
            call do1(s% eLambda_bar)
            if (failed('eLambda_bar')) exit
            call do1(s% rho_bar)
            if (failed('rho_bar')) exit
            call do1(s% Kcond)
            if (failed('Kcond')) exit
            call do1(s% lnK)
            if (failed('lnK')) exit
            call do1(s% dlnK_dlnT)
            if (failed('dlnK_dlnT')) exit
            return
        end do
        ierr = -1
        
    contains
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
        function failed(str)
            character(len=*), intent(in) :: str
            logical :: failed
            failed = .FALSE.
            if (ierr == 0) return
            write (*,*) 'NScool_info_arrays failed for ',trim(str)
            failed = .TRUE.
        end function failed
    end subroutine NScool_info_arrays

    subroutine NScool_work_arrays(s,action,ierr)
        type(NScool_info), pointer :: s
        integer, intent(in) :: action
        integer, intent(out) :: ierr
        integer :: ntab, nz
        
        ierr = 0
        nz = s% nz
        ntab = s% n_tab
        do
            call do1(s% tab_lnT)
            if (failed('tab_lnT')) exit
            call do2(s% lnenu_tab)
            if (failed('lnenu_tab')) exit
            call do2(s% lnCp_tab)
            if (failed('lnCp_tab')) exit
            call do2(s% lnK_tab)
            if (failed('lnK_tab')) exit
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
        function failed(str)
            character(len=*), intent(in) :: str
            logical :: failed
            failed = .FALSE.
            if (ierr == 0) return
            write (*,*) 'NScool_work_arrays failed for ',trim(str)
            failed = .TRUE.
        end function failed
    end subroutine NScool_work_arrays
    
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
