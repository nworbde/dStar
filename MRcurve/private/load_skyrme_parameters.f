module load_skyrme_parameters
	use dStar_core_def
	
contains
	subroutine load_skyrme_table(model,ierr)
		use, intrinsic :: iso_fortran_env, only: error_unit
		character(len=*), intent(in) :: model
		integer, intent(out) :: ierr
		integer :: funit
		type(skyrme), pointer :: tab
		character(len=256) :: model_file
		character(len=32) :: pset
		real(dp) :: gamma, a, b, c, d, e, E16, L, K, E10
		
		model_file = trim(core_datadir)//'/'//trim(model)//'.data'
		open (newunit=funit,file=trim(model_file),status='old',action='read',iostat=ierr)
		if (ierr /= 0) then
			write(error_unit,*) 'unable to open '//trim(model_file)//' for reading'
			return
		end if
		
		! skip the first line
		read(funit,*)
		! read in lines
		call read_one
		if (ierr /= 0) return
		call read_one
		if (ierr /= 0) return
		call read_one
		if (ierr /= 0) return

	contains
		subroutine read_one()
			read (funit,*) pset,gamma,a,b,c,d,e,E16,L,K,E10
			select case (pset)
				case ('matter')
				tab => skyrme_eos(skyrme_matter)
				tab% model = skyrme_matter
				case ('neutron')
				tab => skyrme_eos(skyrme_neutron)
				tab% model = skyrme_neutron
				case ('sym')
				tab => skyrme_eos(skyrme_sym)
				tab% model = skyrme_sym
				case default
				ierr = -1
				write(error_unit,*) 'unable to parse line in '//trim(model)//' for model '//trim(pset)
				return
			end select
			tab% gamma = gamma
			tab% a = a
			tab% b = b
			tab% c = c
			tab% d = d
			tab% e = e
		end subroutine read_one
	end subroutine load_skyrme_table

end module load_skyrme_parameters
