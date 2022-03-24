program run_dStar
	use exceptions_lib
    use NScool_def
    use NScool_lib
    use constants_def, only : boltzmann
    use alt_micro
    use argparse
    
    character(len=*), parameter :: default_inlist_file = 'inlist'
    character(len=64) :: my_dStar_dir
    character(len=64) :: inlist
    real(dp) :: eV_to_MK
    type(NScool_info), pointer :: s=>null()
    integer :: ierr, NScool_id, i
    real(dp), dimension(8) :: pred_Teff, obs_Teff, obs_Teff_sig
    real(dp) :: chi2
    type(failure) :: check_okay=failure(scope='run_dStar')
    
    ierr = 0
    call command_arg_set( &
        & 'dStar_directory',"sets the main dStar root directory",ierr, &
        & flag='D',takes_parameter=.TRUE.)
    if (check_okay% failure(ierr,'set command argument dStar_directory')) stop
    
    call command_arg_set( &
        & 'inlist_file','sets the namelist parameter file',ierr, &
        & flag='I',takes_parameter=.TRUE.)
    if (check_okay% failure(ierr,'set command argument inlist file')) stop
    
    call parse_arguments(ierr)
    if (check_okay% failure(ierr,'parse_arguments')) stop

    my_dStar_dir = trim(get_command_arg('dStar_directory'))
    inlist = trim(get_command_arg('inlist_file'))
    if (len_trim(inlist)==0) inlist = default_inlist_file
 
    call NScool_init(my_dStar_dir, ierr)
    if (check_okay% failure(ierr,'NScool_init')) stop
    
    NScool_id = alloc_NScool(ierr)
    if (check_okay% failure(ierr,'NScool_id')) stop
    
    call NScool_setup(NScool_id,inlist,ierr)
    if (check_okay% failure(ierr,'NScool_setup')) stop
    
    call get_NScool_info_ptr(NScool_id,s,ierr)
    if (check_okay% failure(ierr,'get_NScool_info_ptr')) stop
 
    s% other_set_Qimp => alt_Qimp
    s% other_sf_get_results => alt_sf

    call NScool_create_model(NScool_id,ierr)
    if (check_okay% failure(ierr,'NScool_create_model')) stop

    call NScool_evolve_model(NScool_id,ierr)        
    if (check_okay% failure(ierr,'NScool_evolve_model')) stop
   
    ! we don't want to compare the effective temp. at t = 0, the end of the 
    ! outburst
    pred_Teff = s% Teff_monitor(2:)/1.0e6
    eV_to_MK = 1.602176565e-12_dp/boltzmann/1.0e6
    ! observed effective temperatures (eV) and uncertainties
    obs_Teff = [103.2,88.9,75.5,73.3,71.0,66.0,70.3,63.1] * eV_to_MK
    obs_Teff_sig = [1.7,1.3,2.2,2.3,1.8,4.5,2.1,2.1] * eV_to_MK
    chi2 = sum((pred_Teff-obs_Teff)**2 / obs_Teff_sig**2 )
    
    write(output_unit,*)
    write(output_unit,'(a7,a6,tr3,a12)') 'time','Teff','obs. range'
    write(output_unit,'(a7,a6,tr3,a12)') '[d]','[MK]','[MK]'
    write(output_unit,'(13("-"),tr3,12("-"))')
    do i = 1, 8
        write(output_unit,'(f7.1,f6.3,tr3,f6.3,f6.3)')  &
        & s% t_monitor(i+1), &
        & pred_Teff(i),(obs_Teff(i)-obs_Teff_sig(i)), &
        & (obs_Teff(i)+obs_Teff_sig(i))
    end do
    write(output_unit,*)
    write(output_unit,'(a,f6.2)') 'chi2 = ',chi2
    
    call NScool_shutdown
    
end program run_dStar
