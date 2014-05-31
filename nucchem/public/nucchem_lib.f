module nucchem_lib
    use nucchem_def
    implicit none

    character(len=*), parameter, private :: mod_name = 'nucchem_lib'
    character(len=*), parameter, private :: err_fmt = '(a,":",a,"> ",a)'

contains    
	subroutine nucchem_init(datadir, ierr)
        use iso_fortran_env, only: error_unit
        use nucchem_io
        use nucchem_storage
		character(len=*), intent(in) :: datadir
		integer, intent(out) :: ierr
        character(len=*),parameter :: nuclib_db = 'nuclib_db'
        character(len=160) :: nuclib_filename
        character(len=160) :: nuclib_cache
        integer :: indx

        nuclib_filename = trim(datadir)//'/'//nuclib_db
        nuclib_cache = trim(datadir)//'/cache/'//nuclib_db//'.bin'
        
        ierr = 0
        write(error_unit,'(a)')  &
        & 'loading nuclib from '//trim(nuclib_filename)
        call do_load_nuclib(nuclib_filename,nuclib_cache,ierr)
        write(error_unit,'(/,a,i0,a)')  &
        & 'done. ',nuclib% Nnuclides, &
        & ' nuclides retrieved. now writing nuclide dictionary...'
        call do_parse_nuclides(ierr)
        write(error_unit,'(/,a)') 'done.'

		nucchem_is_initialized = .TRUE.
	end subroutine nucchem_init

	subroutine nucchem_shutdown()
		use utils_lib
		call integer_dict_free(nuclide_dict)
		nucchem_is_initialized = .FALSE.
	end subroutine nucchem_shutdown

	subroutine compute_composition_moments(nnuclides, nuclide_ids, abunds, &
	&   comp, Xsum, ncharged, charged_ids, Yion, abunds_are_mass_fractions, &
	&   exclude_neutrons, renormalize_mass_fractions)
		use iso_fortran_env, only : error_unit
		integer, intent(in) :: nnuclides
		integer, dimension(nnuclides), intent(in) :: nuclide_ids
		real(dp), dimension(nnuclides), intent(in) :: abunds
		type(composition_info_type), intent(out) :: comp
		real(dp), intent(out) :: Xsum
		integer, intent(out) :: ncharged
		integer, dimension(nnuclides), intent(out) :: charged_ids
			! charged_ids(1:ncharged) contain the ids of the charged species
		real(dp), dimension(nnuclides), intent(out) :: Yion
			! Yion(1:ncharged) contains the abundances of the charged species. 
            ! If exclude_neutrons == .TRUE.,
			! then these are renormalized to Y/(1-Yn).
        
            ! optional parameters: defaults are all false
		logical, optional, intent(in) :: abunds_are_mass_fractions
		logical, optional, intent(in) :: exclude_neutrons
			! if true, neutrons will be excluded from the computation of Z, Z53, ...., and Yion will be renorm'd.
		logical, optional, intent(in) :: renormalize_mass_fractions
			! if true, will renormalize the mass fractions

		real(dp), parameter :: fivethird = 5.0_dp/3.0_dp
        real(dp), parameter :: seventhird = 7.0_dp/3.0_dp
		character(len=*), parameter :: routine_name = 'compute_composition_moments'
		real(dp), dimension(nnuclides) :: Z, A, Y
		real(dp), dimension(nnuclides) :: Zion, Aion
		logical :: abnds, excl, rnrm
		logical, dimension(nnuclides) :: ion_mask		
		
        ! set the options
		abnds = .FALSE.; excl = .FALSE.; rnrm = .FALSE.
		if (present(abunds_are_mass_fractions)) abnds = abunds_are_mass_fractions
		if (present(exclude_neutrons)) excl = exclude_neutrons
		if (present(renormalize_mass_fractions)) rnrm = renormalize_mass_fractions
		
		call clear_composition(comp)
		if (.not. nucchem_is_initialized) then
			write(error_unit,err_fmt) mod_name,routine_name, &
			&   'nucchem module is not initialized'
			return
		end if
		if (abnds) then
			Xsum = sum(abunds)
			Y = abunds(1:nnuclides)/nuclib% A(nuclide_ids(1:nnuclides))
		else
			Xsum = dot_product(abunds(1:nnuclides),nuclib% A(nuclide_ids(1:nnuclides)))
			Y = abunds(1:nnuclides)
		end if
		
		if (rnrm) Y = Y/Xsum

		Z = real(nuclib% Z(nuclide_ids(1:nnuclides)),dp)
		A = real(nuclib% A(nuclide_ids(1:nnuclides)),dp)
        
		comp% Ye = dot_product(Y, Z)	! electron fraction is ALWAYS (no. electrons)/(no. baryons)
		comp% Yn = sum(Y, mask = (Z==0.0 .and. A==1.0))

		! mark the charged species
		ion_mask = ( Z > 0.0 .and. A > 0.0 )
		ncharged = count(ion_mask)
		Zion(1:ncharged) = pack(Z, ion_mask)
		Aion(1:ncharged) = pack(A, ion_mask)
		Yion(1:ncharged) = pack(Y, ion_mask)
		charged_ids(1:ncharged) = pack(nuclide_ids, ion_mask)
		
		if (excl) then
			if (comp% Yn > 0.0 .and. comp% Yn /= 1.0) Yion(1:ncharged) = Yion(1:ncharged)/(1.0-comp% Yn)
			call do_compute_moments(Zion(1:ncharged), Aion(1:ncharged), Yion(1:ncharged))
		else
			call do_compute_moments(Z, A, Y)
		end if
				
		contains
		subroutine do_compute_moments(Z,A,Y)
			real(dp), dimension(:), intent(in) :: Z,A,Y
			comp% A = 1.0/sum(min(1.0,max(Y,1.0e-50)))
			comp% Z = dot_product(Y, Z)*comp% A
			comp% Z53 = dot_product(Y, Z**fivethird)*comp% A
			comp% Z2 = dot_product(Y, Z**2)*comp% A
			comp% Z73 = dot_product(Y, Z**seventhird)*comp% A
			comp% Z52 = dot_product(Y, Z**2.5)*comp% A
			comp% ZZ1_32 = dot_product(Y, Z*(Z+1.0)**1.5)*comp% A
			comp% Z2XoA2 = sum(Z**2*Y/A)
			comp% Q = comp% Z2 - (comp% Z)**2
		end subroutine do_compute_moments
	end subroutine compute_composition_moments

	! routines to lookup nuclides
	! 	
	function get_nuclide_index(nuclei) result(indx)
		use iso_fortran_env, only : error_unit
      	use utils_lib, only: integer_dict_lookup
      	character(len=*), intent(in) :: nuclei
      	integer :: indx, ierr
		character(len=*), parameter :: routine_name='get_nuclide_index'
        if (.not. nucchem_is_initialized) then
            write(error_unit,err_fmt) mod_name,routine_name,'nucchem module is not initialized'
            indx = nuclide_not_found
            return
         end if
         ierr = 0
         call integer_dict_lookup(nuclide_dict, nuclei, indx, ierr)
         if (ierr /= 0) indx = nuclide_not_found
  end function get_nuclide_index
  
	! get index of nuclei given charge, neutron numbers
	function get_nuclide_index_from_ZN(Z,N) result(indx)
		use iso_fortran_env, only : error_unit
        use nucchem_def
        
        integer, intent(in) :: Z, N
		character(len=*), parameter :: routine_name='get_nuclide_index'
        character(len=iso_name_length) :: nuclide
        integer :: indx
        
        if (Z < 0 .or. Z > max_element_Z) then
            write(error_unit,err_fmt) mod_name,routine_name,'bad charge number'
            indx = nuclide_not_found
            return
        end if
        if (Z == 0 .and. N == 1) then
            nuclide = 'n'
        else if (Z == 1 .and. N == 0) then
            nuclide = 'p'
        else
            write(nuclide,'(a,i0)') trim(adjustl(element_name(Z))),N+Z
        end if
        
        indx = get_nuclide_index(nuclide)
	end function get_nuclide_index_from_ZN
    
    subroutine clear_composition(comp)
        type(composition_info_type), intent(inout) :: comp
        comp% A = 0.0; comp% Z = 0.0
        comp% Z53 = 0.0; comp% Z2 = 0.0; comp% Z73 = 0.0; comp% Z52 = 0.0
        comp% ZZ1_32 = 0.0; comp% Z2XoA2 = 0.0
        comp% Ye = 0.0; comp% Yn = 0.0; comp% Q = 0.0
    end subroutine clear_composition
    

end module nucchem_lib
