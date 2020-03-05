module PPP_electron
    use constants_def, only: dp
    use conductivity_def

    implicit none
    logical, parameter :: dbg = .TRUE.
    integer, parameter :: NZ = 15, Nrho = 64, NT = 19
    logical :: conductivity_is_initialized = .FALSE.
    integer, parameter :: conductivity_filename_length = 256
    character(len=conductivity_filename_length) :: conductivity_datadir    
    
    character(len=*), parameter :: tablename = 'condall06'
    
    type electron_conductivity_tbl
        logical :: is_loaded
        real(dp) :: lgZmin
        real(dp) :: lgZmax
        real(dp) :: lgrhomin
        real(dp) :: lgrhomax
        real(dp) :: lgTmin
        real(dp) :: lgTmax
        real(dp), dimension(NZ) :: lgZs
        real(dp), dimension(Nrho) :: lgrhos
        real(dp), dimension(NT) :: lgTs
        real(dp), dimension(NT,Nrho,NZ) :: lgK
    end type electron_conductivity_tbl

    type(electron_conductivity_tbl), target :: PPP_tbl

contains
    
    subroutine load_PPP_electron_table(datadir,ierr)
        implicit none
        character(len=*), intent(in) :: datadir
        integer, intent(out) :: ierr
        character(len=*), parameter :: this_routine = 'load_PPP_electron_table'
        character(len=256) :: data_filename, cache_filename
        logical :: have_cache
        type(electron_conductivity_tbl), pointer :: tab
        
        if (dbg) print *,this_routine
        
        ierr = 0
        tab => PPP_tbl
        if (tab% is_loaded) then
            call warning(this_routine//': table is already loaded')
            return
        end if
        
        cache_filename = &
        &   trim(datadir)//'/cache/'//trim(tablename)//'.bin'
        inquire(file=cache_filename,exist=have_cache)
        if (have_cache) then
            call read_PPP_electron_table_cache(cache_filename,tab,ierr)
            if (ierr == 0) return
        end if
        
        ! If we don't have the table cached, or cannot load it, then read the 
        ! datafile and write to cache
        ierr = 0
        data_filename =  &
        &   trim(datadir)//'/'//trim(tablename)//'.dat'
        call read_PPP_electron_table(data_filename,tab,ierr)
        if (failure(this_routine//': unable to load table',ierr)) return
        call write_PPP_electron_table_cache(cache_filename,tab,ierr)
        if (failure(this_routine//': unable to write cache',ierr)) return
    end subroutine load_PPP_electron_table

    subroutine read_PPP_electron_table(datafile,tab,ierr)
        character(len=*), intent(in) :: datafile
        type(electron_conductivity_tbl), pointer :: tab
        integer, intent(out) :: ierr
        character(len=*), parameter :: this_routine='read_PPP_electron_table'
        integer :: unitno, iZ, irho
        real(dp), dimension(NZ) :: Zs
        
        if (dbg) print *,this_routine
        
        open(newunit=unitno,file=datafile,status='old',action='read', &
        &   iostat = ierr)
        if (failure(this_routine//' opening '//datafile,ierr)) then
            close(unitno)
            return
        end if
        
        read(unitno,*)  ! skip first line
        do iZ = 1, NZ
            read(unitno,*,iostat=ierr) Zs(iZ),tab% lgTs(:)
            if (failure(this_routine//': unable to read Zs, Ts',ierr)) then
                close(unitno)
                return
            end if
            do irho = 1, Nrho
                read(unitno,*,iostat=ierr)  &
                &   tab% lgrhos(irho), tab% lgK(:,irho,iZ)
                if (failure(this_routine//': unable to read rho, K',ierr)) then
                    close(unitno)
                    return
                end if
            end do
        end do
        close(unitno)
        tab% lgZs = log10(Zs)
        
        ! set up table values
        tab% lgZmin = minval(tab% lgZs)
        tab% lgZmax = maxval(tab% lgZs)
        tab% lgrhomin = minval(tab% lgrhos)
        tab% lgrhomax = maxval(tab% lgrhos)
        tab% lgTmin = minval(tab% lgTs)
        tab% lgTmax = maxval(tab% lgTs)

        tab% is_loaded = .TRUE.
    end subroutine read_PPP_electron_table
    
    subroutine read_PPP_electron_table_cache(cache_filename,tab,ierr)
        character(len=*), intent(in) :: cache_filename
        type(electron_conductivity_tbl), pointer :: tab
        integer, intent(out) :: ierr
        character(len=*), parameter :: &
        &   this_routine='read_PPP_electron_table_cache'
        integer :: unitno
        
        if (dbg) print *,this_routine
        
        open(newunit=unitno,file=trim(cache_filename),action='read', &
        &   form='unformatted',iostat=ierr)
        if (failure(this_routine//': unable to open '//cache_filename, ierr)) &
        &   return
        
        read(unitno) tab% lgZs
        read(unitno) tab% lgrhos
        read(unitno) tab% lgTs
        read(unitno) tab% lgK
        close(unitno)
        
        ! set up table values
        tab% lgZmin = minval(tab% lgZs)
        tab% lgZmax = maxval(tab% lgZs)
        tab% lgrhomin = minval(tab% lgrhos)
        tab% lgrhomax = maxval(tab% lgrhos)
        tab% lgTmin = minval(tab% lgTs)
        tab% lgTmax = maxval(tab% lgTs)

        tab% is_loaded = .TRUE.
    end subroutine read_PPP_electron_table_cache

    subroutine write_PPP_electron_table_cache(cache_filename,tab,ierr)
        character(len=*), intent(in) :: cache_filename
        type(electron_conductivity_tbl), pointer :: tab
        integer, intent(out) :: ierr
        character(len=*), parameter :: this_routine = &
        &   'write_PPP_electron_table_cache'        
        integer :: unitno
        
        if (dbg) print *,this_routine
        
        open(newunit=unitno,file=trim(cache_filename),action='write', &
        &   form='unformatted',iostat=ierr)
        if (failure(this_routine//': unable to open '//cache_filename,ierr)) &
        &   return
        
        write(unitno) tab% lgZs
        write(unitno) tab% lgrhos
        write(unitno) tab% lgTs
        write(unitno) tab% lgK
        close(unitno)
    end subroutine write_PPP_electron_table_cache

    function failure(msg,ierr)
        use iso_fortran_env, only: error_unit
        character(len=*), intent(in) :: msg
        integer, intent(in) :: ierr
        logical :: failure

        failure = (ierr /= 0)
        if (failure) then
            write(error_unit,'(a)') trim(msg)
        end if
    end function failure
    
    subroutine warning(msg)
        use iso_fortran_env, only: error_unit
        character(len=*), intent(in) :: msg
        write(error_unit,'(a)') trim(msg)
    end subroutine warning
end module PPP_electron
