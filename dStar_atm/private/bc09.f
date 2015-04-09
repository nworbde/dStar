module bc09
    
    ! pointers for the photosphere routine
    integer, parameter :: igrav = 1
    integer, parameter :: itau = 2
    integer, parameter :: iTeff = 3
    integer, parameter :: number_photosphere_rpar = 3
    integer, parameter :: ihandle = 1
    integer, parameter :: number_photosphere_ipar = 1

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
        real(dp), intent(in) :: Teff,grav,tau
        real(dp), intent(out) :: Pphoto
        integer, intent(in) :: eos_handle
        integer, intent(out) :: ierr
        
    end subroutine find_photospheric_pressure    

    real(dp) function photosphere(rho, dfdrho, lrpar, rpar, lipar, ipar, ierr)
       ! returns with ierr = 0 if was able to evaluate f and df/dx at x
       ! if df/dx not available, it is okay to set it to 0
       use const_def, only: dp
       integer, intent(in) :: lrpar, lipar
       real(dp), intent(in) :: rho
       real(dp), intent(out) :: dfdrho
       integer, intent(inout), pointer :: ipar(:) ! (lipar)
       real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
       integer, intent(out) :: ierr
       real(dp) :: gravity, tau_ph, Teff, P, kappa
       integer :: eos_handle
       type(composition_info_type) :: ion
       type(conductivity_components) :: K
       real(dp), dimension(num_dStar_eos_results) :: res
       integer :: phase, ncharged
       integer, dimension(:) :: charged_ids
       
       ierr = 0
       gravity = rpar(igrav)
       tau_ph = rpar(itau)
       Teff = rpar(iTeff)
       eos_handle = ipar(ihandle)
       
       call eval_crust_eos(eos_handle,rho,Teff,ionic,ncharged,charged_ids,Yion, &
       		&   res,phase,use_default_nuclear_size)
       
       P = exp(res(i_lnP))

       get_thermal_conductivity(rho,T,chi,Gamma,eta,ionic,K,which_components)
       kappa = K(icond_kap)
       
       return P - 2.0*onethird*gravity/kappa
    end function f


end module bc09
