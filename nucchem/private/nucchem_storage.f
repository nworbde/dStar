module nucchem_storage

contains
    subroutine allocate_nuclib_data(w,n,ierr)
        use nucchem_def, only: nuclib_data, npfcn
        type(nuclib_data), intent(out) :: w
        integer, intent(in) :: n
        integer, intent(out) :: ierr
        
        ierr = 0
        w% Nnuclides = 0
        allocate(w% name(n), w% provenance(n), w% A(n), w% Z(n), w% N(n),  &
        & w% spin(n), w% mass_excess(n), w% pfcn(npfcn,n), stat=ierr)
        if (ierr /= 0) return
        w% Nnuclides = n
    end subroutine allocate_nuclib_data
    
    subroutine free_nuclib_data(w)
        use nucchem_def, only: nuclib_data
        type(nuclib_data), intent(inout) :: w
        if (allocated(w% name)) deallocate(w% name)
        if (allocated(w% provenance)) deallocate(w% provenance)
        if (allocated(w% A)) deallocate(w% A)
        if (allocated(w% Z)) deallocate(w% Z)
        if (allocated(w% N)) deallocate(w% N)
        if (allocated(w% spin)) deallocate(w% spin)
        if (allocated(w% mass_excess)) deallocate(w% mass_excess)
        if (allocated(w% pfcn)) deallocate(w% pfcn)
        w% Nnuclides = 0
    end subroutine free_nuclib_data
    

end module nucchem_storage
