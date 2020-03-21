module PPP_electron
    use constants_def, only: dp
    use conductivity_def

    implicit none
    ! error codes
    integer, parameter :: unable_to_load_table = -1
    integer, parameter :: unable_to_compute_interpolation = -2
    integer, parameter :: unable_to_write_cache = -3
    integer, parameter :: unable_to_evaluate_conductivity = -4
    integer, parameter :: table_is_not_loaded = 1
    integer, parameter, private :: NZ = 15, Nrho = 64, NT = 19
    integer, parameter, private :: conductivity_filename_length = 256
    character(len=conductivity_filename_length) :: conductivity_datadir    
    
    character(len=*), parameter :: tablename = 'condall06'
    
    type electron_conductivity_tbl
        logical :: is_loaded = .FALSE.
        integer :: linear_T
        integer :: linear_rho
        real(dp) :: Zmin
        real(dp) :: Zmax
        real(dp) :: rhomin
        real(dp) :: rhomax
        real(dp) :: Tmin
        real(dp) :: Tmax
        real(dp), dimension(NZ) :: lgZs
        real(dp), dimension(Nrho) :: lgrhos
        real(dp), dimension(NT) :: lgTs
        real(dp), dimension(NT,Nrho,NZ) :: lgK    ! table and coeff.
    end type electron_conductivity_tbl

    type(electron_conductivity_tbl), target, save :: PPP_tbl
    ! can't put this in the table because of the target attribute
    real(dp), dimension(4,NT,Nrho,NZ), target, private, save :: workspace

contains
    
    subroutine load_PPP_electron_table(datadir,ierr)
        use exceptions_lib
        implicit none
        character(len=*), intent(in) :: datadir
        integer, intent(out) :: ierr
        type(electron_conductivity_tbl), pointer :: tab
        character(len=*), parameter :: this_routine = 'load_PPP_electron_table'
        character(len=256) :: data_filename, cache_filename
        logical :: have_cache
        type(alert) :: status=alert(scope=this_routine)
        type(failure) :: read_table_error=failure(scope=this_routine, &
        &   message='unable to load table')
        type(failure) :: write_cache_error=failure(scope=this_routine, &
        &   message='unable to write cache')
        
        ierr = 0
        tab => PPP_tbl
        if (tab% is_loaded) then
            call status% report('table is already loaded')
            return
        end if
        
        cache_filename = trim(datadir)//'/cache/'//trim(tablename)//'.bin'
        inquire(file=cache_filename,exist=have_cache)
        if (have_cache) then
            call read_PPP_electron_table_cache(cache_filename,tab,ierr)
            if (ierr == 0) then
                tab% is_loaded = .TRUE.
                return
            end if
        end if
        
        ! If we don't have the table cached, or cannot load it, then read the 
        ! datafile and write to cache
        ierr = 0
        data_filename = trim(datadir)//'/'//trim(tablename)//'.dat'
        call read_PPP_electron_table(data_filename,tab,ierr)
        if (read_table_error% raised(ierr)) then
            ierr = unable_to_load_table
            return
        end if
        tab% is_loaded = .TRUE.

        call write_PPP_electron_table_cache(cache_filename,tab,ierr)
        if (write_cache_error% raised(ierr)) then
            ierr = unable_to_write_cache
            return
        end if
    end subroutine load_PPP_electron_table
    
    subroutine free_PPP_electron_table()
        PPP_tbl% is_loaded = .FALSE.
    end subroutine free_PPP_electron_table

    subroutine construct_interpolation_coefficients(ierr)
        use interp_2d_lib_db, only: interp_mkbicub_db
        use exceptions_lib
        implicit none
        integer, intent(out) :: ierr
        real(dp), dimension(:), pointer :: ftab=>null()
        type(electron_conductivity_tbl), pointer :: tab
        character(len=*), parameter :: &
        &    this_routine='construct_interpolation_coefficients'
        integer, parameter :: not_a_knot = 0
        integer :: iZ
        real(dp), dimension(NT) :: bcrhomin, bcrhomax
        real(dp), dimension(Nrho) :: bcTmin, bcTmax
        type(failure) :: interpolation_error=failure(scope=this_routine, &
        &   message='unable to construct interpolation')
        
        tab => PPP_tbl
        workspace(1,:,:,:) = tab% lgK
        bcrhomin = 0.0_dp; bcrhomax= 0.0_dp
        bcTmin = 0.0_dp; bcTmax = 0.0_dp
        
        do iZ = 1, NZ
            ftab(1:4*Nrho*NT) => workspace(:,:,:,iZ)
            call interp_mkbicub_db( &
            &   tab% lgTs,NT,tab% lgrhos,Nrho,ftab,NT, &
            &   not_a_knot,bcTmin,not_a_knot,bcTmax, &
            &   not_a_knot,bcrhomin,not_a_knot,bcrhomax, &
            &   tab% linear_T,tab% linear_rho,ierr)
            
            if (interpolation_error% raised(ierr)) then
                ierr = unable_to_compute_interpolation
                tab% is_loaded = .FALSE.
                return
            end if
        end do

    end subroutine construct_interpolation_coefficients

    subroutine read_PPP_electron_table(datafile,tab,ierr)
        use exceptions_lib
        implicit none
        character(len=*), intent(in) :: datafile
        type(electron_conductivity_tbl), pointer :: tab
        integer, intent(out) :: ierr
        character(len=*), parameter :: this_routine='read_PPP_electron_table'
        integer :: unitno, iZ, irho
        real(dp), dimension(NZ) :: Zs
        type(failure) :: io_error=failure(scope=this_routine)
        
        open(newunit=unitno,file=datafile,status='old',action='read', &
        &   iostat = ierr)
        if (io_error% raised(ierr,'opening '//trim(datafile))) then
            close(unitno)
            return
        end if
        
        read(unitno,*)  ! skip first line
        do iZ = 1, NZ
            read(unitno,*,iostat=ierr) Zs(iZ),tab% lgTs(:)
            if (io_error% raised(ierr,'unable to read Zs, Ts')) then
                close(unitno)
                return
            end if
            do irho = 1, Nrho
                read(unitno,*,iostat=ierr)  &
                &   tab% lgrhos(irho), tab% lgK(:,irho,iZ)
                if (io_error% raised(ierr,'unable to read rho, K')) then
                    close(unitno)
                    return
                end if
            end do
        end do
        close(unitno)
        
        tab% lgZs = log10(Zs)
        call set_table_limits(tab)
    end subroutine read_PPP_electron_table
    
    subroutine read_PPP_electron_table_cache(cache_filename,tab,ierr)
        use exceptions_lib
        implicit none
        character(len=*), intent(in) :: cache_filename
        type(electron_conductivity_tbl), pointer :: tab
        integer, intent(out) :: ierr
        character(len=*), parameter :: &
        &   this_routine='read_PPP_electron_table_cache'
        integer :: unitno
        type(failure) :: io_error=failure(scope=this_routine)
        
        open(newunit=unitno,file=trim(cache_filename),action='read', &
        &   form='unformatted',iostat=ierr)
        if (io_error% raised(ierr,'unable to open '//cache_filename)) return
        
        read(unitno) tab% lgZs
        read(unitno) tab% lgTs
        read(unitno) tab% lgrhos
        read(unitno) tab% lgK
        call set_table_limits(tab)
        close(unitno)
    end subroutine read_PPP_electron_table_cache

    subroutine write_PPP_electron_table_cache(cache_filename,tab,ierr)
        use exceptions_lib
        implicit none
        character(len=*), intent(in) :: cache_filename
        type(electron_conductivity_tbl), pointer :: tab
        integer, intent(out) :: ierr
        character(len=*), parameter :: this_routine = &
        &   'write_PPP_electron_table_cache'        
        integer :: unitno
        type(failure) :: io_error=failure(scope=this_routine)
        
        open(newunit=unitno,file=trim(cache_filename),action='write', &
        &   form='unformatted',iostat=ierr)
        if (io_error% raised(ierr,'unable to open '//cache_filename)) return
        
        write(unitno) tab% lgZs
        write(unitno) tab% lgTs
        write(unitno) tab% lgrhos
        write(unitno) tab% lgK
        close(unitno)
    end subroutine write_PPP_electron_table_cache
    
    subroutine set_table_limits(tab)
        implicit none
        type(electron_conductivity_tbl), pointer :: tab
        
        ! set up table markers
        tab% Zmin = 10.0_dp**minval(tab% lgZs)
        tab% Zmax = 10.0_dp**maxval(tab% lgZs)
        tab% Tmin = 10.0_dp**(minval(tab% lgTs))
        tab% Tmax = 10.0_dp**(maxval(tab% lgTs))
        tab% rhomin = 10.0_dp**(minval(tab% lgrhos))
        tab% rhomax = 10.0_dp**(maxval(tab% lgrhos))        
    end subroutine set_table_limits
    
    subroutine eval_PPP_electron_table(rho,T,Z,K,ierr)
        use iso_fortran_env, only: error_unit
        use num_lib, only: binary_search
        use interp_2d_lib_db, only: interp_evbicub_db
        use exceptions_lib
        implicit none
        real(dp), intent(in) :: rho,T,Z
        real(dp), intent(out) :: K
        integer, intent(out) :: ierr
        character(len=*), parameter :: this_routine='eval_PPP_electron_table'
        type(electron_conductivity_tbl), pointer :: tab
        real(dp) :: lgrho, lgT, lgZ, lgK0, lgK1, lgZ0, lgZ1
        real(dp), dimension(:), pointer :: ftab=>null()
        integer, dimension(6) :: ict
        real(dp), dimension(6) :: lgK
        integer :: lZ
        type(warning) :: loaded_table=warning(scope=this_routine, &
        &   message='table is not loaded')
        type(failure) :: interpolation_error=failure(scope=this_routine)
        
        ierr = 0
        K = 0.0_dp
        tab => PPP_tbl
        if (.not.tab% is_loaded) then
            ierr = table_is_not_loaded
            if (loaded_table% raised(ierr)) return
        end if
            
        lgrho = log10(clip_to_table(rho, tab% rhomin, tab% rhomax))
        lgT = log10(clip_to_table(T, tab% Tmin, tab% Tmax))
        lgZ = log10(clip_to_table(Z, tab% Zmin, tab% Zmax))
        ict = [ 1, 0, 0, 0, 0, 0]

        ! search in Z
        lZ = binary_search(NZ,tab% lgZs,-1,lgZ)
        lgZ0 = tab% lgZs(lZ)
        ftab(1:4*Nrho*NT) => workspace(:,:,:,lZ)
        call interp_evbicub_db(lgT,lgrho,tab% lgTs,NT,tab% lgrhos,Nrho, &
        &   tab% linear_T,tab% linear_rho,ftab,NT,ict,lgK,ierr)
        if (interpolation_error% raised(ierr,'interpolation Z0')) then
            ierr = unable_to_evaluate_conductivity
            return
        end if
        lgK0 = lgK(1)
        ! catch the edge case
        if (lZ == NZ) then
            K = 10.0_dp**lgK0
            return
        end if
        lZ = lZ + 1
        lgZ1 = tab% lgZs(lZ)
        ftab(1:4*Nrho*NT) => workspace(:,:,:,lZ)
        call interp_evbicub_db(lgT,lgrho,tab% lgTs,NT,tab% lgrhos,Nrho, &
        &   tab% linear_T,tab% linear_rho,ftab,NT,ict,lgK,ierr)
        if (interpolation_error% raised(ierr,'interpolation Z1')) then
            ierr = unable_to_evaluate_conductivity
            return
        end if
        lgK1 = lgK(1)
        
        ! linear interpolate in Z
        K = 10.0_dp**(lgK0*(lgZ1-lgZ)/(lgZ1-lgZ0)+lgK1*(lgZ-lgZ0)/(lgZ1-lgZ0))
        
    contains
        function clip_to_table(x,xmin,xmax) result(xc)
            real(dp), intent(in) :: x, xmin, xmax
            real(dp) :: xc
            xc = min(max(x,xmin),xmax)
        end function clip_to_table
    end subroutine eval_PPP_electron_table

end module PPP_electron
