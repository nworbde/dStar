program generate_HZ90_table
    use const_def, only: dp
    use exceptions_lib
    use nucchem_def
    use nucchem_lib
    use composition_handler

    integer, parameter :: HZ90_number_table_points = 2048
    real(dp), parameter :: HZ90_lgPmin = 22.0
    real(dp), parameter :: HZ90_lgPmax = 33.5
    real(dp), parameter :: HZ90_temperature = 1.0e8_dp
    character(len=*), parameter :: datadir = '../data/'
    character(len=*), parameter :: cache_stem = 'HZ90'
    
    integer, parameter :: HZ90_number_isos = 19
    character(len=iso_name_length), parameter, dimension(HZ90_number_isos) :: HZ90_network = &
    &   [ character(len=iso_name_length) :: &
    &   'n', &
    &   'mg40', &
    &   'mg44', &
    &   'mg48', &
    &   'si46', &
    &   'si50', &
    &   'si54', &
    &   's52', &
    &   's56', &
    &   's60', &
    &   'ar56', &
    &   'ar62', &
    &   'ar66', &
    &   'ca56', &
    &   'ca68', &
    &   'ti56', &
    &   'ti88',  &
    &   'cr56', &
    &   'fe56' ]
    integer, parameter :: HZ90_number_transitions = 17
    integer, parameter :: HZ90_number_layers = HZ90_number_transitions+1
    character(len=iso_name_length), parameter, dimension(HZ90_number_layers) :: HZ90_ion_composition = [ &
    &   character(len=iso_name_length) :: &
    &   'fe56','cr56','ti56','ca56','ar56','s52','si46','mg40','ca68','ar62','s56','si50', &
    &   'mg44','ar66','s60','si54','mg48','ti88' ]
    real(dp),parameter, dimension(HZ90_number_transitions) :: HZ90_transition_pressures = [ &
    &   7.235d26,9.569d27,1.152d29,4.747d29,1.361d30,1.980d30,2.253d30, &
    &   2.637d30,2.771d30,3.216d30,3.825d30,4.699d30,6.043d30,7.233d30,9.238d30, &
    &   1.228d31,1.602d31 ]
    real(dp), parameter, dimension(HZ90_number_layers) :: Xn = [ &
    &   0.0, 0.0, 0.0, 0.0, 0.0,  &
    &   4.0/56.0, 10.0/56.0, 16.0/56.0, 22.0/56.0, &
    &   50.0/112.0, 56.0/112.0, 62.0/112.0, 68.0/112.0, 158.0/224.0, &
    &   164.0/224.0, 170.0/224.0, 176.0/224.0, 360.0/448.0 ]
    
    real(dp), dimension(HZ90_number_table_points) :: lgP
    character(len=iso_name_length), dimension(HZ90_number_isos) :: isos
    real(dp), dimension(HZ90_number_isos,HZ90_number_table_points) :: Y
    
    integer, dimension(max_nnuclib) :: network_indcs
    integer :: nz, nion, i, indx, n_indx, ierr
    integer, dimension(HZ90_number_isos) :: indcs
    real(dp), dimension(HZ90_number_transitions) :: lgPt
    real(dp) :: delta_lgP
    type(assertion) :: nucchem_is_online = assertion(scope='generate_HZ90_table', &
    & message='unable to initialize nucchem module')
    type(assertion) :: write_cache=assertion(scope='generate_HZ90_table', &
    & message='unable to write cache')
    type(alert) :: status=alert(scope='generate_HZ90_table')
    character(len=128) :: alert_msg
    character(len=64) :: cache_filename

    ierr = 0
    call nucchem_init('../../data/',ierr)
    call nucchem_is_online% assert(ierr==0)
    
    ! local variables
    lgPt = log10(HZ90_transition_pressures)
    nz = HZ90_number_table_points
    nion = HZ90_number_isos
    
    ! indcs are locations of nuclei in nuclib database
    indcs = [(get_nuclide_index(HZ90_network(i)),i=1,HZ90_number_isos)]
    ! network_indcs allow for reverse lookup of specific nuclei
    network_indcs = 0
    network_indcs(indcs) = [(i,i=1,HZ90_number_isos)]
    n_indx = network_indcs(get_nuclide_index('n'))
    
    ! allocate the pressure table
    delta_lgP = HZ90_lgPmax-HZ90_lgPmin
    lgP = [ (HZ90_lgPmin + real(i-1,dp)*(delta_lgP)/real(HZ90_number_table_points-1,dp),  &
        & i = 1,HZ90_number_table_points)]
    
    ! loop over layers to set composition
    Y = 0.0
    ! first layer
    indx = network_indcs(get_nuclide_index(HZ90_ion_composition(1)))
    where(lgP <= lgPt(1)) 
        Y(n_indx,:) = Xn(1)
        Y(indx,:) = (1.0-Xn(1))/nuclib% A(indcs(indx))
    end where
    ! middle layers
    do i = 2, HZ90_number_transitions
        indx = network_indcs(get_nuclide_index(HZ90_ion_composition(i)))
        where(lgP > lgPt(i-1) .and. lgP <= lgPt(i))
            Y(n_indx,:) = Xn(i)
            Y(indx,:) = (1.0-Xn(i))/nuclib% A(indcs(indx))
        end where
    end do
    ! last layer
    indx = network_indcs(get_nuclide_index(HZ90_ion_composition(HZ90_number_layers)))
    where (lgP > lgPt(HZ90_number_transitions))
        Y(n_indx,:) = Xn(HZ90_number_layers)
        Y(indx,:) = (1.0-Xn(HZ90_number_layers))/nuclib% A(indcs(indx))
    end where

    cache_filename = datadir//cache_stem//'.bin'
    call write_composition_cache(cache_filename,nz,nion,HZ90_network,HZ90_temperature,lgP,Y,ierr)
    call write_cache% assert(ierr==0)
    write(alert_msg,'(a,a)') 'wrote HZ90 table to ',cache_filename
    call status% report(alert_msg)
end program generate_HZ90_table

