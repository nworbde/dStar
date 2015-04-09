module bc09

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

end module bc09
