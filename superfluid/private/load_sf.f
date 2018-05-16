module load_sf
    use superfluid_def
    
    contains
    subroutine load_sf_table(prefix,sf_type,ierr)
        use, intrinsic :: iso_fortran_env, only: error_unit
        use interp_1d_def
        use interp_1d_lib
        
        character(len=*), intent(in) :: prefix
        integer, intent(in) :: sf_type
        integer, intent(out) :: ierr
        character(len=256) :: filename
        integer :: funit
        type(sf_table_type), pointer :: tab
        real(dp), dimension(:), pointer :: work=>null()
        real(dp), pointer, dimension(:,:) :: fval
        
        if (sf_type <= max_number_sf_types .and. sf_type > 0) then
            tab => sf_tables(sf_type)
        else
            ierr = -9
            write (error_unit,*) 'load_sf_table: malformed superfluid type'
            return
        end if
        
        ! if the tables are already allocated, issue a warning and scrub the table
        if (tab% is_loaded) then
            write(error_unit,*) &
            &    'load_sf_table WARNING: overwriting already loaded table'
            call free_sf_table(tab)
        end if
        
        call get_file_name(prefix,sf_type,filename,ierr)
        
        open(newunit=funit,file=filename,status='old',action='read', &
        & iostat = ierr)
        if (ierr /= 0) then
            write (error_unit,*) 'unable to open file for reading'
            return
        end if
        
        ! skip the first 3 lines
        read(funit,*)
        read(funit,*)
        read(funit,*)
        
        ! now get the number of points and the minimum, maximum kF
        read(funit,*) tab% nv, tab% kF_min, tab% kF_max
        allocate(tab% kF(tab% nv),tab% f(4*tab% nv),stat=ierr)
        if (ierr /= 0 ) then
            write(error_unit,*) 'load_sf_table: unable to allocate memory'
            close(funit)
            return
        end if
        
        fval(1:4,1:tab% nv) => tab% f(1:4*tab% nv)
        read(funit,*) tab% kF(:)
        read(funit,*) fval(1,:)
        
        close(funit)
        
        allocate(work(tab% nv*pm_work_size))
        call interp_pm(tab% kF, tab% nv, tab% f, pm_work_size, work, &
        &    'load_sf_table', ierr)
        deallocate(work)
        tab% is_loaded = .TRUE.
    end subroutine load_sf_table

    subroutine get_file_name(prefix,sf_type,filename,ierr)
        use, intrinsic :: iso_fortran_env, only: error_unit
        character(len=*), intent(in) :: prefix
        integer, intent(in) :: sf_type
        character(len=*), intent(out) :: filename
        integer, intent(out) :: ierr
        character(len=2) :: typecode
        
        ierr = 0
        filename = ''
        select case (sf_type)
            case (proton_1S0)
                typecode = 'p1'
            case (neutron_1S0)
                typecode = 'n1'
            case (neutron_3P2)
                typecode = 'n3'
            case default
                ierr = -1
                write(error_unit,*) 'get_file_name: invalid superfluid type'
        end select
        
        if (ierr == 0) then
            write (filename,'(a)') trim(sf_datadir) // '/' //trim(prefix) //  &
                &   '_' // typecode // '.data'
        end if
    end subroutine get_file_name

    subroutine free_sf_table(tab)
        type(sf_table_type), pointer :: tab
        tab% which_gap = undefined_gap
        tab% nv = 0
        tab% kF_min = 0.0
        tab% kF_max = 0.0
        if (allocated(tab% kF)) deallocate(tab% kF)
        if (associated(tab% f)) then
            deallocate(tab% f)
            nullify(tab% f)
        end if
        tab% is_loaded = .FALSE.
    end subroutine free_sf_table

end module load_sf
