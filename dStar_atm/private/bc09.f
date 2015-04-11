module bc09
    
    ! pointers for the photosphere routine
    integer, parameter :: igrav = 1
    integer, parameter :: itau = 2
    integer, parameter :: iTeff = 3
    integer, parameter :: iPph = 4
    integer, parameter :: iKph = 5
    integer, parameter :: number_photosphere_rpar = 5
    integer, parameter :: ihandle = 1
    integer, parameter :: ichem_id = 2
    integer, parameter :: number_photosphere_ipar = 2

contains
    
	subroutine do_get_bc09_Teff(grav, Plight, Tb, Teff, flux)
		use constants_def
		real(dp), intent(in) :: grav	! surface gravity, in the local frame
		real(dp), intent(in) :: Plight	! pressure at which layer of light elements terminates
		real(dp), intent(in), dimension(:) :: Tb	! temperature at a base column
		real(dp), intent(out), dimension(:) :: Teff, flux	! effective temperature and flux
		real(dp) :: eta, g14
		real(dp), dimension(size(Tb)) :: Tb9, Teff6_4
        
        ! make a very dense table of Tb(Teff); then interpolate to get Teff(Tb)
        integer :: size_tab ! = 4.0*size(Tb)
        real(dp), dimension(:), allocatable :: tabTb9, tabTeff_4
        
        ! compute dense table
        
        ! interpolate from dense table to get finished product
    end subroutine do_get_bc09_Teff
    
    subroutine find_photospheric_pressure(Teff,grav,tau,Pphoto,eos_handle,ierr)
        use constants_def
        use nucchem_def, only : nuclide_not_found
    	use nucchem_lib, only : get_nuclide_index
        
        real(dp), intent(in) :: Teff,grav,tau
        real(dp), intent(out) :: Pphoto
        integer, intent(in) :: eos_handle
        integer, intent(out) :: ierr
        real(dp) :: root_ph,rho,dfdrho
        integer, pointer :: ipar(:) => null() ! (lipar)
        real(dp), pointer :: rpar(:) => null()  ! (lrpar)
        integer :: i
        
        ierr = 0
        allocate(ipar(number_photosphere_ipar), rpar(number_photosphere_rpar))

        ! hardwire photosphere compostion to pure He
        ipar(ichem_id) = get_nuclide_index('he4')
        if (ipar(ichem_id) == nuclide_not_found) then
            ierr = nuclide_not_found
            return
        end if
                
        rpar(igrav) = grav
        rpar(iTeff) = Teff
        rpar(itau) = tau
        ipar(ihandle) = eos_handle
        
        do i = 1, 20
            rho = 10.0**(-1+2.0*(i-1)/19.0)
            root_ph = photosphere(rho,dfdrho,number_photosphere_rpar,rpar, &
                & number_photosphere_ipar,ipar,ierr)
            print '(2(a,es11.4))','rho = ',rho,'P - 2g/3k = ',root_ph   
            Pphoto = rpar(iPph)
            print '(a,es11.4)','Pphoto = ',Pphoto
        end do
        deallocate(ipar, rpar)
    end subroutine find_photospheric_pressure    

    real(dp) function photosphere(rho, dfdrho, lrpar, rpar, lipar, ipar, ierr)
       ! returns with ierr = 0 if was able to evaluate f and df/dx at x
       ! if df/dx not available, it is okay to set it to 0
       use constants_def
       use nucchem_def
       use nucchem_lib
       use dStar_eos_lib
       use conductivity_lib
       
       integer, intent(in) :: lrpar, lipar
       real(dp), intent(in) :: rho
       real(dp), intent(out) :: dfdrho
       integer, intent(inout), pointer :: ipar(:) ! (lipar)
       real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
       integer, intent(out) :: ierr
       real(dp) :: gravity, tau_ph, Teff, P, kappa, Gamma, eta, Xsum, chi
       integer :: eos_handle
       type(composition_info_type) :: ionic
       type(conductivity_components) :: K
       real(dp), dimension(num_dStar_eos_results) :: res
       integer :: phase, ncharged
       integer, dimension(1) :: charged_ids, chem_ids
	   real(dp),dimension(1) :: Y,Yion
	   
       ierr = 0
       
       ! check inputs
       if (rho < 0.0) then
           ierr = -9
           return
       end if
       
       gravity = rpar(igrav)
       tau_ph = rpar(itau)
       Teff = rpar(iTeff)
       eos_handle = ipar(ihandle)
	   chem_ids = ipar(ichem_id)
       chi = use_default_nuclear_size
       ! single species only
       Y(1) = 1.0/nuclib% A(chem_ids(1))

       call compute_composition_moments(1,chem_ids,Y,ionic,Xsum,ncharged, charged_ids, Yion,  &
   			& exclude_neutrons = .TRUE.)
       call eval_crust_eos(eos_handle,rho,Teff,ionic,ncharged,charged_ids,Yion, &
       		&   res,phase,chi)
       
       P = exp(res(i_lnP))
       rpar(iPph) = P
       eta = res(i_Theta) !1.0/TpT
       Gamma = res(i_Gamma)
       call get_thermal_conductivity(rho,Teff,chi, &
           & Gamma,eta,ionic,K,which_components=cond_use_only_kap)
       kappa = 4.0*onethird*arad*clight*Teff**3/rho/K% kap
       rpar(iKph) = kappa
	   dfdrho = 0.0
       
       photosphere = P - 2.0*onethird*gravity/kappa
    end function photosphere

end module bc09
