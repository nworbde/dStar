module bc09

contains
    
    subroutine find_photospheric_pressure(Teff,grav,tau,Pphoto,eos_handle,ierr)
        use constants_def
        real(dp), intent(in) :: Teff,grav,tau
        real(dp), intent(out) :: Pphoto
        integer, intent(in) :: eos_handle
        integer, intent(out) :: ierr
        
        
    end subroutine find_photospheric_pressure    

end module bc09
