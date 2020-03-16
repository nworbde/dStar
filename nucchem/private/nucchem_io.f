module nucchem_io

    ! the winvne format
    character(len=*), parameter :: pfcn_fmt = '(8(f9.2))'
    character(len=*), parameter :: iso_fmt = '(a5)'
    character(len=*), parameter :: line0_fmt = '(a5,f12.3,2i4,f6.1,f10.3,a6)'

contains
    
    subroutine do_load_nuclib(filename,cache_filename,ierr)
        use, intrinsic :: iso_fortran_env, only: iostat_end, error_unit
        use exceptions_lib
        use nucchem_def
        use nucchem_storage
        
        character(len=*), intent(in) :: filename,cache_filename
        integer, intent(out) :: ierr
        type(nuclib_data) :: tmp_n
        integer :: i, iso_count,ios, nuclib_unitno
        real, dimension(24) :: pfcn
        logical :: have_cache
        type(failure) :: load_failure=failure(scope='do_load_nuclib')
        
        ierr = 0
        inquire(file=cache_filename,exist=have_cache)
        if (have_cache) then
            call do_read_nuclib_cache(cache_filename,ierr)
            if (ierr == 0) return
        end if
        
        open(newunit=nuclib_unitno, file=trim(filename), iostat=ierr, &
        & status="old", action="read")
        if (load_failure% raised(ierr,message='opening '//trim(filename))) &
        &   return

        ! allocate a temporary to hold the library
        call allocate_nuclib_data(tmp_n,max_nnuclib,ierr)
        if (load_failure% raised(ierr,message='allocating storage')) return

        i = 1
        do
            read(nuclib_unitno,line0_fmt,iostat=ios) &
            & tmp_n% name(i), tmp_n% A(i), tmp_n% Z(i), tmp_n% N(i),  &
            & tmp_n% spin(i), tmp_n% mass_excess(i),tmp_n% provenance(i)
            if (ios == iostat_end) exit
            read(nuclib_unitno,pfcn_fmt) tmp_n% pfcn(1:8,i)
            read(nuclib_unitno,pfcn_fmt) tmp_n% pfcn(9:16,i)
            read(nuclib_unitno,pfcn_fmt) tmp_n% pfcn(17:24,i)
            tmp_n% name(i) = adjustl(tmp_n% name(i))
            i = i+1
            if (mod(i,100) == 0)  &
            & write (error_unit,'(a)',advance='no') '.'
        end do
        close(nuclib_unitno)

        tmp_n% Nnuclides = i-1

        call copy_nuclib_data(tmp_n,ierr)
        if (load_failure% raised(ierr,'copying nuclib data')) return
        call free_nuclib_data(tmp_n)
        
        if (.not.have_cache) then
            call do_write_nuclib_cache(cache_filename,ierr)
        end if
    contains
        
        subroutine copy_nuclib_data(old_data,ierr)
            type(nuclib_data), intent(in) :: old_data
            integer, intent(out) :: ierr
            integer :: n
        
            if (nuclib% Nnuclides /= 0) call free_nuclib_data(nuclib)
            call allocate_nuclib_data(nuclib,old_data% Nnuclides, ierr)
            if (ierr /= 0) return
        
            n = old_data% Nnuclides
            nuclib% name(:) = old_data% name(1:n)
            nuclib% provenance(:) = old_data% provenance(1:n)
            nuclib% A(:) = old_data% A(1:n)
            nuclib% Z(:) = old_data% Z(1:n)
            nuclib% N(:) = old_data% N(1:n)
            nuclib% spin(:) = old_data% spin(1:n)
            nuclib% mass_excess(:) = old_data% mass_excess(1:n)
            nuclib% pfcn(:,:) = old_data% pfcn(:,1:n)
        end subroutine copy_nuclib_data
        
    end subroutine do_load_nuclib    

    subroutine do_read_nuclib_cache(filename,ierr)
        use exceptions_lib
        use nucchem_def
        use nucchem_storage
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ierr
        integer :: cache_unitno, n
        type(failure) :: cache_failure=failure(scope='do_read_nuclib_cache')
        
        ierr = 0
        
        open(newunit=cache_unitno, file=trim(filename), iostat=ierr, &
        &   action="read",status="old",form="unformatted")
        if (cache_failure% raised(ierr,message='opening '//trim(filename))) &
        &    return
        
        read(cache_unitno) n
        call allocate_nuclib_data(nuclib,n,ierr)
        if (cache_failure% raised(ierr,message='allocating nuclib')) return
        read(cache_unitno) nuclib% name
        read(cache_unitno) nuclib% provenance
        read(cache_unitno) nuclib% A
        read(cache_unitno) nuclib% Z
        read(cache_unitno) nuclib% N
        read(cache_unitno) nuclib% spin
        read(cache_unitno) nuclib% mass_excess
        read(cache_unitno) nuclib% pfcn
        read(cache_unitno) pfcn_T9
        close(cache_unitno)
    end subroutine do_read_nuclib_cache
     
    subroutine do_write_nuclib_cache(filename,ierr)
        use exceptions_lib
        use nucchem_def
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ierr
        integer :: cache_unitno
        type(failure) :: cache_failure=failure(scope='do_write_nuclib_cache')
        
        ierr = 0
        if (nuclib% Nnuclides == 0) return
        
        open(newunit=cache_unitno, file=trim(filename), iostat=ierr, &
        &   action="write",form="unformatted")
        if (cache_failure% raised(ierr, &
        &   message='opening '//trim(filename))) return
        
        write(cache_unitno) nuclib% Nnuclides
        write(cache_unitno) nuclib% name
        write(cache_unitno) nuclib% provenance
        write(cache_unitno) nuclib% A
        write(cache_unitno) nuclib% Z
        write(cache_unitno) nuclib% N
        write(cache_unitno) nuclib% spin
        write(cache_unitno) nuclib% mass_excess
        write(cache_unitno) nuclib% pfcn
        write(cache_unitno) pfcn_T9

        close(cache_unitno)
    end subroutine do_write_nuclib_cache
   
    subroutine do_parse_nuclides(ierr)
        use utils_def, only: integer_dict
        use utils_lib, only: integer_dict_define, &
        & integer_dict_lookup, integer_dict_free, integer_dict_create_hash
        use nucchem_def
        
        integer, intent(out) :: ierr
        integer :: indx, ikey
         
        if (associated(nuclide_dict)) call integer_dict_free(nuclide_dict)
        do indx = 1, nuclib% Nnuclides
            call integer_dict_define(nuclide_dict, &
            & trim(adjustl(nuclib% name(indx))),indx,ierr)
        end do
        call integer_dict_create_hash(nuclide_dict,ierr)
    end subroutine do_parse_nuclides
    
!     function io_failure(action,ierr)
!         use iso_fortran_env, only: error_unit
!         character(len=*), intent(in) :: action
!         integer, intent(in) :: ierr
!         logical :: io_failure
!         character(len=*), parameter :: err_format = &
!         &  '("Error while ",a,". error code = ",i0)'
!
!         if (ierr == 0) then
!             io_failure = .FALSE.
!             return
!         end if
!         io_failure = .TRUE.
!         write (error_unit,err_format) action, ierr
!     end function io_failure

end module nucchem_io
