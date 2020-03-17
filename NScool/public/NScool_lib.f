module NScool_lib

    logical, save :: dbg = .FALSE.

contains
    subroutine NScool_init(my_dStar_dir,ierr)
        use init, only: do_NScool_init
        character(len=*), intent(in) :: my_dStar_dir
        integer, intent(out) :: ierr
        ierr = 0
        call do_NScool_init(my_dStar_dir, ierr)
    end subroutine NScool_init
    
    subroutine NScool_shutdown()
       use init, only: do_NScool_shutdown
        call do_NScool_shutdown
    end subroutine NScool_shutdown
    
    function alloc_NScool(ierr)
        use init, only: alloc_NScool_data
        integer, intent(out) :: ierr
        integer :: alloc_NScool
        alloc_NScool = alloc_NScool_data(ierr)
    end function alloc_NScool
    
    subroutine dealloc_NScool(id,ierr)
        use init, only: dealloc_NScool_data
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        call dealloc_NScool_data(id,ierr)
    end subroutine dealloc_NScool

    ! read a namelist and set parameters
    subroutine NScool_setup(id,inlist,ierr)
       use NScool_ctrls_io, only: do_one_setup
       integer, intent(in) :: id
       character(len=*), intent(in) :: inlist
       integer, intent(out) :: ierr
       call do_one_setup(id,inlist,ierr)
    end subroutine NScool_setup
   
    subroutine NScool_create_model(id,ierr)
       use init, only: model_builder
       integer, intent(in) :: id
       integer, intent(out) :: ierr

       call model_builder(id, ierr)
    end subroutine NScool_create_model

    subroutine NScool_evolve_model(id, ierr)
        use constants_def, only : julian_day, amu, ergs_to_mev
        use exceptions_lib
        use NScool_def, only : NScool_info, get_NScool_info_ptr
        use NScool_evolve, only: do_integrate_crust
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        type(NScool_info), pointer :: s
        integer :: i_epoch
        type(assertion) :: got_pointer=assertion(scope='NScool_evolve_model', &
        &   message='unable to get pointer')
        type(failure) :: integration_error=failure(scope='NScool_evolve_model', &
        &   message='integration error; exiting loop over epochs')

        ierr = 0
        call get_NScool_info_ptr(id,s,ierr)
        call got_pointer% assert(ierr == 0)

        do i_epoch = 1, s% number_epochs
            s% epoch_start_time = s% epoch_boundaries(i_epoch-1)
            s% epoch_duration = s% epoch_boundaries(i_epoch) - s% epoch_start_time
            s% epoch_id = i_epoch
            s% Mdot = s% epoch_Mdots(i_epoch)
            if (i_epoch > 1) then
                if (s% suppress_first_step_output) then
                    s% starting_number_for_profile = s% model
                else
                    s% starting_number_for_profile = s% model + 1
                end if
            end if
            call do_integrate_crust(id,ierr)
            if (integration_error% raised(ierr)) exit
            s% t_monitor(i_epoch) = s% tsec / julian_day
            s% Teff_monitor(i_epoch) = s% Teff * s% ePhi(1)
            ! need to check about the frame in which Mdot is defined
            if (s% Mdot > 0.0) then
                s% Qb_monitor(i_epoch) = s% Lsurf/s% Mdot * amu*ergs_to_mev
            else
                s% Qb_monitor(i_epoch) = 0.0
            end if
        end do
    end subroutine NScool_evolve_model

end module NScool_lib
