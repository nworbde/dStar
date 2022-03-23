program process_abuntime
    use exceptions_lib
    use math_lib
    use constants_def
    use constants_lib
    use nucchem_def
    use nucchem_lib
	use superfluid_def
	use superfluid_lib
	use dStar_eos_lib
    use abuntime
    use composition_handler

    character(len=*), parameter :: cache_dir = '../data/'
    integer :: argument_count
    character(len=64) :: argument
    character(len=64) :: abuntime_stem,abuntime_filename, abuntime_cache
    real(dp) :: lgP_increment, abundance_threshold, T
    integer :: nz, nion !, ncharged
    real(dp), dimension(:,:), allocatable :: Yion,Yout,Ytot
    real(dp), dimension(:), allocatable :: lgP,lgPtot
    character(len=iso_name_length), dimension(:), allocatable :: isos,isonet, &
    &   network
    integer :: ierr
    integer :: nnet, i, j, ntot,nztot
    
    character(len=128) :: alert_msg
    type(alert) :: status=alert(scope='process_abuntime')
    type(assertion) :: command_arguments=assertion(scope='process_abuntime', &
        & message='USAGE: process_abuntime <file stem> <increment in lgP> <abundance threshold>')
    type(assertion) :: init_okay=assertion(scope='process_abuntime', &
        & message='unable to initialize constants')
    type(assertion) :: nucchem_load_okay=assertion(scope='process_abuntime', &
        & message='unable to initialize nucchem')
    type(assertion) :: abuntime_load_okay=assertion(scope='process_abuntime', &
        & message='read_abuntime file')
    type(assertion) :: abuntime_reduction_okay=assertion(scope='process_abuntime', &
        & message='reduce abuntime file')
    type(assertion) :: abuntime_extension_okay=assertion(scope='process_abuntime', &
        & message='extend abuntime file')
    type(assertion) :: abuntime_cache_okay=assertion(scope='process_abuntime', &
        & message='write abuntime cache')

    argument_count = command_argument_count()
    call command_arguments% assert(argument_count==3)
    
    call get_command_argument(1,argument)
    read(argument,*) abuntime_stem
    call get_command_argument(2,argument)
    read(argument,*) lgP_increment
    call get_command_argument(3,argument)
    read(argument,*) abundance_threshold
    
    abuntime_filename = cache_dir//trim(abuntime_stem)
    abuntime_cache = cache_dir//trim(abuntime_stem)//'.bin'

    ierr = 0
    call math_init()
    call constants_init('../..','',ierr)
    call init_okay% assert(ierr==0)
    
    call nucchem_init(ierr)
    call nucchem_load_okay% assert(ierr==0)
    
    write(alert_msg,'(a)') 'reading '//trim(abuntime_filename)
    call status% report(alert_msg)
    call read_abuntime(trim(abuntime_filename), nz, nion, T, lgP, &
    &   isos, Yion, ierr,lgP_increment)
    call abuntime_load_okay% assert(ierr==0)

    write(alert_msg,'(a)') 'processing table...'
    call reduce_abuntime(nz,nion,isos,lgP,Yion,nnet,isonet,Yout,ierr,abundance_threshold)
    call abuntime_reduction_okay% assert(ierr == 0)
    
    call extend_abuntime(nz,nnet,isonet,lgP,Yout,lgP_increment,nztot,ntot,network,lgPtot,Ytot,ierr)
    call abuntime_extension_okay% assert(ierr==0)
    
    write(alert_msg,'(a)') 'writing cache to '//trim(abuntime_cache)
    call status% report(alert_msg)
    call write_composition_cache(trim(abuntime_cache),nztot,ntot,network,lgPtot,Ytot,ierr)
    call abuntime_cache_okay% assert(ierr == 0)

    deallocate(lgP,isos,Yion,isonet,Yout)
	call nucchem_shutdown

end program process_abuntime
