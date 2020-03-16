module dStar_eos_private_def
    use dStar_eos_def
    
    ! error codes
    character(len=80), dimension(-3:0), parameter :: dstar_eos_private_def_errors =  &
    &   [character(len=80):: &
    &   'invalid dStar_eos_handle', &
    &   'broken handle for dStar_eos_id', &
    &   'no available crust eos handle', &
    &   'no errors']

    type dStar_eos_general_info
        real(dp) :: Gamma_melt
        real(dp) :: rsi_melt
        real(dp) :: Ythresh
        real(dp) :: pasta_transition
        real(dp) :: cluster_transition
        logical :: suppress_warnings
        integer :: handle
        logical :: in_use
    end type dStar_eos_general_info
    
    integer, parameter :: max_dStar_eos_handles = 10
    type(dStar_eos_general_info), dimension(max_dStar_eos_handles), target ::  dStar_eos_handles

contains
    subroutine dStar_eos_def_init()
        integer :: i
        do i = 1, max_dStar_eos_handles
            dStar_eos_handles(i)% Gamma_melt = default_Gamma_melt
            dStar_eos_handles(i)% rsi_melt = default_rsi_melt
            dStar_eos_handles(i)% Ythresh = default_Y_threshold
            dStar_eos_handles(i)% pasta_transition = default_pasta_transition
            dStar_eos_handles(i)% cluster_transition = default_cluster_transition
            dStar_eos_handles(i)% suppress_warnings = default_suppress_eos_warnings
            dStar_eos_handles(i)% handle = i
            dStar_eos_handles(i)% in_use = .FALSE.
        end do
    end subroutine dStar_eos_def_init

    function do_alloc_dStar_eos(ierr) result(dStar_eos_id)
        integer, intent(out) :: ierr
        integer :: dStar_eos_id, i
        ierr = 0
        dStar_eos_id = -1
        do i = 1, max_dStar_eos_handles
            if (.not. dStar_eos_handles(i)% in_use) then
                dStar_eos_handles(i)% in_use = .TRUE.
                dStar_eos_id = i
                exit
            end if
        end do
        
        if (dStar_eos_id == -1) then
            ierr = -1
            return
        end if
        if (dStar_eos_handles(i)% handle /= dStar_eos_id) then
            ierr = -2
            return
        end if
    end function do_alloc_dStar_eos
    
    subroutine do_free_dStar_eos(handle)
        integer, intent(in) :: handle
        if (handle >= 1 .and. handle <= max_dStar_eos_handles)  &
        &   dStar_eos_handles(handle)% in_use = .FALSE.
    end subroutine do_free_dStar_eos
    
    subroutine get_dStar_eos_ptr(handle,rq,ierr)
        integer, intent(in) :: handle
        type(dStar_eos_general_info), pointer, intent(out) :: rq
        integer, intent(out) :: ierr
        if (handle < 1 .or. handle > max_dStar_eos_handles) then
            ierr = -1
            rq => null()
            return
        end if
        rq => dStar_eos_handles(handle)
        ierr = 0
    end subroutine get_dStar_eos_ptr
    
end module dStar_eos_private_def
