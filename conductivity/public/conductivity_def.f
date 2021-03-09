module conductivity_def
    use constants_def, only: dp
    
    ! flag to control which ee scattering fmla. is used. 
    ! default: Shternin & Yakovlev (2006)
    integer, parameter :: icond_sy06 = 1
    ! Potekhin, Chabrier and Yakovlev (1997)
    integer, parameter :: icond_pcy = 2

    ! flag to control which eQ scattering fmla. is used.
    ! default: Potekhin (private comm.)
    integer, parameter :: icond_eQ_potekhin = 1
    ! Dany Page (private comm.)
    integer, parameter :: icond_eQ_page = 2

    type conductivity_components
        real(dp) :: total
        real(dp) :: ee
        real(dp) :: ei
        real(dp) :: eQ
        real(dp) :: sf
        real(dp) :: nQ
        real(dp) :: np
        real(dp) :: kap
        real(dp) :: electron_total
        real(dp) :: neutron_total
    end type conductivity_components
    
    type conductivity_general_info
        ! flags to turn off some of the more exploratory transport mechanisms
        logical :: include_neutrons
        logical :: include_superfluid_phonons
        integer :: ee_scattering_fmla
        integer :: eQ_scattering_fmla
        ! radiative opacity ramps up linearly between full_off and full_on
        real(dp) :: rad_full_off_lgrho
        real(dp) :: rad_full_on_lgrho
        ! defines region where we linearly blend from the table to analytical treatment
        real(dp) :: tab_off_lgrho
        real(dp) :: tab_on_lgrho
        integer :: handle
        logical :: in_use
    end type conductivity_general_info
    
    integer, parameter :: max_conductivity_handles = 10
    type (conductivity_general_info), target :: &
    &   conductivity_handles(max_conductivity_handles)
    
    real(dp), parameter :: default_rad_full_off_lgrho = 10.0_dp
    real(dp), parameter :: default_rad_full_on_lgrho = 9.0_dp
    real(dp), parameter :: default_tab_on_lgrho = 9.0_dp
    real(dp), parameter :: default_tab_off_lgrho = 9.0_dp
    
    logical, save :: conductivity_is_initialized = .FALSE.
    
contains
    
    subroutine conductivity_def_init()
        integer :: i
        do i = 1, max_conductivity_handles
            conductivity_handles(i)% include_neutrons = .TRUE.
            conductivity_handles(i)% include_superfluid_phonons = .TRUE.
            conductivity_handles(i)% ee_scattering_fmla = icond_sy06
            conductivity_handles(i)% eQ_scattering_fmla = icond_eQ_potekhin
            conductivity_handles(i)% rad_full_off_lgrho = default_rad_full_off_lgrho
            conductivity_handles(i)% rad_full_on_lgrho = default_rad_full_on_lgrho
            conductivity_handles(i)% tab_off_lgrho = default_tab_off_lgrho
            conductivity_handles(i)% tab_on_lgrho = default_tab_on_lgrho
            conductivity_handles(i)% handle = i
            conductivity_handles(i)% in_use = .FALSE.
        end do    
    end subroutine conductivity_def_init
    
    function do_alloc_conductivity(ierr) result(conductivity_id)
        integer, intent(out) :: ierr
        integer :: conductivity_id, i
        ierr = 0
        conductivity_id = -1
        do i = 1, max_conductivity_handles
            if (.not.conductivity_handles(i)% in_use) then
                conductivity_handles(i)% in_use = .TRUE.
                conductivity_id = i
                exit
            end if
        end do
        
        if (conductivity_id == -1) then
            ierr = -1
            return
        end if
        if (conductivity_handles(i)% handle /= conductivity_id) then
            ierr = -2
            return
        end if
    end function do_alloc_conductivity
    
    subroutine do_free_conductivity(handle)
        integer, intent(in) :: handle
        if (handle >= 1 .or. handle <= max_conductivity_handles) &
        &   conductivity_handles(handle)% in_use = .FALSE.
    end subroutine do_free_conductivity
    
    subroutine get_conductivity_ptr(handle,rq,ierr)
        integer, intent(in) :: handle
        type(conductivity_general_info), pointer, intent(out) :: rq
        integer, intent(out) :: ierr
        if (handle < 1 .or. handle > max_conductivity_handles) then
            ierr = -1
            rq => null()
            return
        end if
        rq => conductivity_handles(handle)
        ierr = 0
    end subroutine get_conductivity_ptr
    
end module conductivity_def
