program run_dStar
	use iso_fortran_env, only : output_unit, error_unit
    use NScool_def
    use NScool_lib
    use constants_def, only : boltzmann
    use alt_micro
    
    character(len=*), parameter :: my_dStar_dir = '/path/to/local/dStar'
    character(len=*), parameter :: inlist = 'inlist'
    real(dp) :: eV_to_MK
    type(NScool_info), pointer :: s=>null()
    integer :: ierr, NScool_id, i
    real(dp), dimension(8) :: pred_Teff, obs_Teff, obs_Teff_sig
    real(dp) :: chi2
    
    ierr = 0
    call NScool_init(my_dStar_dir, ierr)
    call check_okay('NScool_init',ierr)
    
    NScool_id = alloc_NScool(ierr)
    call check_okay('NScool_id',ierr)
    
    call NScool_setup(NScool_id,inlist,ierr)
    call check_okay('NScool_setup',ierr)
    
    call get_NScool_info_ptr(NScool_id,s,ierr)
    call check_okay('get_NScool_info_ptr',ierr)
 
    s% other_set_Qimp => alt_Qimp
    s% other_sf_get_results => alt_sf

    call NScool_create_model(NScool_id,ierr)
    call check_okay('NScool_create_model',ierr)

    call NScool_evolve_model(NScool_id,ierr)        
    call check_okay('NScool_evolve_model',ierr)
   
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
    
contains
	subroutine check_okay(msg,ierr)
		character(len=*), intent(in) :: msg
		integer, intent(inout) :: ierr
		if (ierr /= 0) then
			write (error_unit,*) trim(msg)//': ierr = ',ierr
			if (ierr < 0) stop
		end if
	end subroutine check_okay
end program run_dStar
