module composition_handler

    implicit none
    
contains
    subroutine write_composition_cache(cache_filename,nz,nion,isos,lgP,Y,ierr)
        use const_def, only: dp
        use exceptions_lib
        use nucchem_def, only: iso_name_length
        character(len=*), intent(in) :: cache_filename
        integer, intent(in) :: nz, nion
        character(len=iso_name_length), intent(in), dimension(:) :: isos ! nion
        real(dp), intent(in), dimension(:) :: lgP   ! nz
        real(dp), intent(in), dimension(:,:) :: Y   ! nion, nz
        integer, intent(out) :: ierr
        integer :: unitno
        type(failure) :: open_cache_failure=failure(scope='write_composition_cache')
    
        ierr = 0
        open(newunit=unitno, file=trim(cache_filename),action='write', &
            &   form='unformatted',iostat=ierr)
        if (open_cache_failure% raised(ierr,message='opening '//trim(cache_filename))) return

        write(unitno) nz
        write(unitno) nion
        write(unitno) isos
        write(unitno) lgP
        write(unitno) Y
        close(unitno)
    end subroutine write_composition_cache

    subroutine read_composition_cache(cache_filename,nz,nion,isos,lgP,Y,ierr)
        use const_def, only: dp
        use exceptions_lib
        use nucchem_def, only: iso_name_length
        character(len=*), intent(in) :: cache_filename
        integer, intent(out) :: nz,nion
        character(len=iso_name_length), dimension(:), allocatable, intent(out) :: isos
        real(dp), dimension(:), allocatable, intent(out) :: lgP
        real(dp), dimension(:,:), allocatable, intent(out) :: Y
        integer, intent(out) :: ierr
        integer :: unitno
        type(failure) :: open_cache_failure=failure(scope='read_composition_cache')
        type(failure) :: allocate_table_failure=failure(scope='read_composition_cache', &
            & message='allocating composition tables')

        open(newunit=unitno,file=trim(cache_filename), &
        &   action='read',status='old',form='unformatted', iostat=ierr)
        if (open_cache_failure% raised(ierr,message='opening '//trim(cache_filename))) return

        read(unitno) nz
        read(unitno) nion
        allocate(isos(nion), lgP(nz), Y(nion,nz),stat=ierr)
        if (allocate_table_failure% raised(ierr)) then
            close(unitno)
            return
        end if
        read(unitno) isos
        read(unitno) lgP
        read(unitno) Y
        close(unitno)
    end subroutine read_composition_cache
    
end module composition_handler
