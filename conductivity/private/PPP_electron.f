module PPP_electron
    use constants_def, only: dp
    use conductivity_def

    implicit none
    integer, parameter :: NZ = 15, Nrho = 64, NT = 19
    character(len=*), parameter :: conductivity_file_basename = 'condall06'

    ! error codes
    integer, parameter :: PPP_table_not_loaded = -1
    
    type electron_conductivity_tbl
        logical :: loaded
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
    
    subroutine read_PPP_electron_table(datadir,ierr)
        character(len=*), intent(in) :: datadir
        integer, intent(out) :: ierr
        character(len=256) :: datafile, cachefile
        integer :: iounit
        
        type(electron_conductivity_tbl), pointer :: rq
        
        rq => PPP_tbl
        
        ! Attempt read from cache
        
        ! If unable, read from datafile and write to cache
        
        ! Bail
            
    end subroutine read_PPP_electron_table

    subroutine write_PPP_electron_table_cache(datadir, ierr)
        character(len=*), intent(in) :: datadir
        integer, intent(out) :: ierr
        character(len=*), parameter :: this_routine = &
        &   'write_PPP_electron_table_cache'
        type(electron_conductivity_tbl), pointer :: rq
        integer :: iounit
        character(len=256) :: cache_name
        
        rq => PPP_tbl
        
        if (.not. rq% loaded) then
            ierr = PPP_table_not_loaded
            call failure(ierr,this_routine)
            return
        endif
    end subroutine write_PPP_electron_table_cache

    subroutine failure(ierr,routine)
        use iso_fortran_env, only: error_unit
        integer :: intent(in) :: ierr
        character(len=*), intent(in) :: routine
        select case(ierr)
        case(PPP_table_not_loaded)
            write(error_unit,'(a)') routine//': PPP table is not loaded'
        end select
end module PPP_electron
