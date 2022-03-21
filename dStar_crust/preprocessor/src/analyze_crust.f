program analyze_crust
    use exceptions_lib
    use math_lib
    use constants_def, only: dp
    use constants_lib
    use nucchem_def
    use nucchem_lib
    use superfluid_lib
    use dStar_eos_lib
    use dStar_crust_def
    use dStar_crust_lib
    use composition_handler
    use dStar_crust_mod, only : find_densities
    
    integer :: i, ierr
    real(dp), dimension(:), allocatable :: lgP, lgRho, lgEps
    integer, dimension(:), allocatable :: charged_ids, indcs
    character(len=iso_name_length), dimension(:), allocatable :: isos
    real(dp), dimension(:,:), allocatable :: Y, Yion
    integer :: ncharged, N, Nisos
    character(len=*), parameter :: datadir = '../../data/crust_data/'
    character(len=128) :: composition_filename, output_filename, composition_stem
    type(composition_info_type), dimension(:), allocatable :: ion_info
    integer :: eos_handle, iounit, argument_count
    real(dp) :: Tref, Xsum
    character(len=*), parameter :: routine_name='analyze_crust'
    type(assertion) :: command_arguments=assertion(scope=routine_name, &
        & message='USAGE: '//routine_name//' <file stem>')
    type(failure) :: read_composition_error=failure(scope=routine_name, &
    &   message='unable to read composition')
    type(failure) :: allocation_error=failure(scope=routine_name)
    type(failure) :: interpolation_error=failure(scope=routine_name)
    
    argument_count = command_argument_count()
    call command_arguments% assert(argument_count==1)

    call get_command_argument(1,composition_stem)
    composition_filename = datadir//trim(composition_stem)//'.bin'
    output_filename = trim(composition_stem)//'-table'
    
    call math_init()
    call constants_init('',ierr)
    call check_okay('constants_init',ierr)
    
    call nucchem_init('../../data',ierr)
    call check_okay('nucchem_init',ierr)
    
    call sf_startup('../../data',ierr)
    call check_okay('sf_startup',ierr)
    
    call sf_load_gaps('ns','gc','t72',ierr)
    call check_okay('sf_load_gaps',ierr)
    
    call dStar_eos_startup('../../data')
    call check_okay('dStar_eos_startup',ierr)
    
    eos_handle = alloc_dStar_eos_handle(ierr)
    call check_okay('alloc_dStar_eos_handle',ierr)
    
    ! switch off the warnings about quantum effects
    call dStar_eos_set_controls(eos_handle,suppress_warnings=.TRUE.)
    
    call dStar_crust_startup('../../data',ierr)
    call check_okay('dStar_crust_startup',ierr)

    Tref = 1.0d8

    ierr = 0
    call read_composition_cache(composition_filename,N,Nisos,isos,lgP,Y,ierr)
    if (read_composition_error% raised(ierr)) return
    
    allocate(lgRho(N), lgEps(N), Yion(Nisos,N), ion_info(N), indcs(Nisos), charged_ids(Nisos), stat=ierr)
    if (allocation_error% raised(ierr,'unable to allocate storage')) return
    ! set up markers for the nuclei
    indcs = [(get_nuclide_index(adjustl(isos(i))),i=1,Nisos)]

    ! set the composition moments and densities
    do i = 1, N
        call compute_composition_moments(Nisos, indcs, Y(:,i), &
        &   ion_info(i), Xsum, ncharged, charged_ids, Yion(:,i), &
        &   abunds_are_mass_fractions=.FALSE., exclude_neutrons=.TRUE.)
    end do
    call find_densities(eos_handle,lgP,lgRho,lgEps,Yion,ncharged,charged_ids,ion_info,Tref)
    
    ! output results
    open(newunit=iounit,file=output_filename,action='write',iostat=ierr)
    do i = 1, N
        write(iounit,'(7(f6.3,tr1))') lgP(i), lgEps(i), lgRho(i),  &
        &   ion_info(i)% Z, ion_info(i)% A, ion_info(i)% Ye, ion_info(i)% Yn
    end do
    close(iounit)
    call dStar_crust_shutdown
    call dStar_eos_shutdown
    call sf_shutdown
    call nucchem_shutdown

contains
    subroutine check_okay(msg,ierr)
        character(len=*), intent(in) :: msg
        integer, intent(inout) :: ierr
        type(assertion) :: okay=assertion(scope='main')
        
        call okay% set_message(msg)
        call okay% assert(ierr == 0)
    end subroutine check_okay

end program analyze_crust
