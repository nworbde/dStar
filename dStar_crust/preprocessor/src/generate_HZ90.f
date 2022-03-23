program generate_HZ90_table
    use math_lib
    use const_def, only: dp
    use exceptions_lib
    use constants_lib
    use nucchem_def
    use nucchem_lib
    use HZ90_comp
    use composition_handler

    implicit none
    integer, parameter :: HZ90_number_table_points = 2048
    real(dp), parameter :: HZ90_lgPmin = 22.0_dp
    real(dp), parameter :: HZ90_lgPmax = 33.5_dp
    character(len=*), parameter :: cache_dir = '../data/'
    character(len=*), parameter :: cache_stem = 'HZ90'
    
    real(dp), dimension(HZ90_number_table_points) :: lgP
    character(len=iso_name_length), dimension(HZ90_number_isos) :: isos
    real(dp), dimension(HZ90_number_isos,HZ90_number_table_points) :: Y
    integer :: ierr, nz, nion, i
    real(dp) :: delta_lgP
    type(assertion) :: nucchem_is_online = assertion(scope='generate_HZ90_table', &
    & message='unable to initialize nucchem module')
    type(assertion) :: write_cache=assertion(scope='generate_HZ90_table', &
    & message='unable to write cache')
    type(alert) :: status=alert(scope='generate_HZ90_table')
    character(len=128) :: alert_msg
    character(len=64) :: cache_filename

    ierr = 0
    call math_init()
    call constants_init('../..','',ierr)
    call nucchem_init(ierr)
    call nucchem_is_online% assert(ierr==0)
    
    nz = HZ90_number_table_points
    nion = HZ90_number_isos
    
    ! allocate the pressure table
    delta_lgP = HZ90_lgPmax-HZ90_lgPmin
    lgP = [ (HZ90_lgPmin + real(i-1,dp)*(delta_lgP)/real(HZ90_number_table_points-1,dp),  &
        & i = 1,HZ90_number_table_points)]
    
    call do_generate_HZ90_table(lgP,Y)
    
    cache_filename = cache_dir//cache_stem//'.bin'
    call write_composition_cache(cache_filename,nz,nion,HZ90_network,lgP,Y,ierr)
    call write_cache% assert(ierr==0)
    write(alert_msg,'(a,a)') 'wrote HZ90 table to ',cache_filename
    call status% report(alert_msg)
end program generate_HZ90_table

