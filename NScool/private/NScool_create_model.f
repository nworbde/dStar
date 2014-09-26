module create_model
    use NScool_def
    use constants_def
    use constants_lib
    use nucchem_def
    use nucchem_lib
    use superfluid_lib
    use dStar_eos_lib
    use dStar_crust_def
    use dStar_crust_lib
    use NScool_crust_tov


contains
    
    subroutine do_setup_crust_zones(s, ierr)
        use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
        use storage
        type(NScool_info), pointer :: s
        type(tov_model_type), pointer :: stov
        integer, intent(out) :: ierr
        real(dp), dimension(:), pointer :: y
    
        write (error_unit,*) 'establishing crust zones...'
        
        write (error_unit,*) 'initializing microphysics...'
        call nucchem_init('../../data',ierr)
        if (failure('nucchem_init')) return
	
        call sf_startup('../../data',ierr)
        if (failure('sf_startup')) return
	
        ! gaps are hard-wired for now
        call sf_load_gaps('ns','gc','t72',ierr)
        if (failure('sf_load_gaps')) return
	
        call dStar_eos_startup('../../data')
        if (failure('dStar_eos_startup')) return
	
        s% eos_handle = alloc_dStar_eos_handle(ierr)
        if (failure('alloc_dStar_eos_handle')) return
    
        ! switch off the warnings about quantum effects
        call dStar_eos_set_controls(s% eos_handle,suppress_warnings=.TRUE.)
    
        write (error_unit,*) 'loading crust model...'
        call dStar_crust_startup('../../data',ierr)
        if (failure('dStar_atm_startup')) return
        
        call dStar_crust_load_table('hz90',s% eos_handle, s% Tcore, ierr)
        if (failure('dStar_crust_load_table')) return

        write(error_unit,*) 'integrating TOV equations...'
        allocate(y(num_tov_variables))

        call tov_integrate(log10(s% Pcore), log10(s% Psurf), s% Mcore, s% Rcore, s% target_resolution_lnP, y, ierr)
        if (failure('tov_integrate')) return
        deallocate(y)
        
        ! allocate the info arrays
        stov => tov_model
        s% nz = stov% nzs
        s% nisos = dStar_crust_get_composition_size()
        call allocate_NScool_info_arrays(s, ierr)
        
        ! facial information
        s% P_bar(1:s% nz) = stov% pressure(stov% nzs:1:-1) * pressure_g
        s% ePhi_bar(1:s% nz) = exp(stov% potential(stov% nzs:1:-1))
        s% m(1:s% nz) = stov% baryon(stov% nzs:1:-1) * mass_g
        s% eLambda_bar(1:s% nz) = 1.0/sqrt(1.0-2.0*(stov% mass(stov% nzs:1:-1)+s% Mcore)/ &
        &   (stov% radius(stov% nzs:1:-1) + s% Rcore*1.0e5/length_g))
        
        ! body information
        s% dm(1:s% nz-1) = s% m(1:s% nz-1) - s% m(2:s% nz)
        s% rho(1:s% nz-1) = s% dm(1:s% nz-1)/(stov% volume(stov% nzs:2:-1) - stov% volume(stov% nzs-1:1:-1))/length_g**3
        ! interpolate metric functions
        s% eLambda(1:s% nz-1) = 0.5*(s% eLambda_bar(1:s% nz-1) + s% eLambda_bar(2:s% nz))
        s% ePhi(1:s% nz-1) = 0.5*(s% ePhi_bar(1:s% nz-1) + s% ePhi_bar(2:s% nz))
        
    contains
        function failure(str)
            character(len=*), intent(in) :: str
            logical :: failure
            failure = .FALSE.
            if (ierr == 0) return
            write (error_unit,*) 'do_setup_crust_zones failure in ',trim(str)
            failure = .TRUE.
        end function failure
        
    end subroutine do_setup_crust_zones
    
    subroutine do_setup_crust_composition(s, ierr)
        use storage
        type(NScool_info), pointer :: s
        integer, intent(out) :: ierr
        real(dp), allocatable, dimension(:) :: lgP_bar
        
        ! this routine assumes that we have already run do_setup_crust_zones

        allocate(lgP_bar(s% nz))
        lgP_bar = log10(s% P_bar)
        
        call allocate_NScool_iso_arrays(s, ierr)
        call dStar_crust_get_composition(lgP_bar, s% ncharged, s% charged_ids, s% Yion_bar, s% Xneut_bar, s% ionic_bar, ierr)
        
        ! for the cells, *for now*, inherit composition of the bottom face
        s% Yion(1:s% ncharged, 1:s% nz-1) = s% Yion_bar(1:s% ncharged, 2:s% nz)
        s% ionic(1:s% nz-1) = s% ionic_bar(2:s% nz)
        
        deallocate(lgP_bar)
    end subroutine do_setup_crust_composition
    
    

end module create_model
