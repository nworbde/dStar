module skyrme
	! implements the Skyrme parameterizations from B A Brown.
	!
	use constants_def
	use dStar_eos_def

contains

	subroutine eval_skyrme_eos(id,n,u,p,dur,dpr)
		integer, intent(in) :: id
		real(dp), intent(in) :: n	! density of relevant constituent, fm**-3
		real(dp), intent(out) :: u	! internal energy, MeV per constituent
		real(dp), intent(out) :: p	! pressure, MeV fm**-3
		real(dp), intent(out) :: dur	! du/dlnn
		real(dp), intent(out) :: dpr	! dp/dlnn
		type(skyrme_parameter_set_type), pointer :: s
		
		s => skyrme_eos(id)
        u = s% a*n + s% b*n**s% gamma + s% c*n**twothird +  &
        	& s% d*n**fivethird + s% e*n**3
		dur =  s% a*n + s% b*s% gamma*n**s% gamma +  &
	        	& twothird*s% c*n**twothird + fivethird*s% d*n**fivethird &
	        	& +3.0*s% e*n**3
		p = n*(s% a*n + s% b*s% gamma*n**s% gamma +  &
			& twothird*s% c*n**twothird + fivethird*s% d*n**fivethird &
			& +3.0*s% e*n**3)
		dpr = 2*s% a*n + s% b*s% gamma*(s% gamma+1.0)*n**s% gamma  &
			& + s% c*twothird*fivethird*n**twothird  &
			& + s% d*fivethird*seventhird*n**fivethird + 12.0*s% e*n**3
		dpr = dpr*n
	end subroutine eval_skyrme_eos

	subroutine load_skyrme_table(model,eos_datadir,ierr)
		use, intrinsic :: iso_fortran_env, only: error_unit
		character(len=*), intent(in) :: model,eos_datadir
		integer, intent(out) :: ierr
		integer :: funit
		type(skyrme_parameter_set_type), pointer :: tab
		character(len=256) :: model_file
		character(len=32) :: pset
		character(len=18) :: parameter_set_type
		real(dp) :: gamma, a, b, c, d, e, E16, L, K, E10
		
		model_file = trim(eos_datadir)//'/'//trim(model)//'.data'
		open (newunit=funit,file=trim(model_file),status='old',action='read', &
		&	iostat=ierr)
		if (ierr /= 0) then
			write(error_unit,*)  &
				& 'unable to open '//trim(model_file)//' for reading'
			return
		end if
		
		! skip the first line
		read(funit,*)
		! read in tag for parameter set
		read(funit,'(a18)') parameter_set_type
		! skip next line
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
			a = 0.0; b = 0.0; c = 0.0; d = 0.0; e = 0.0; gamma = 0.0
			select case (parameter_set_type)
				case('4 parameter Skyrme')
				read (funit,*) pset,gamma,a,b,c,d,E16,L,K,E10
				case('5 parameter Skyrme')
				read (funit,*) pset,gamma,a,b,c,d,e,E16,L,K,E10
				case default
				ierr = -6
				write(error_unit,*) 'unrecognized type of parameter set'
				return
			end select
			select case (pset)
				case ('matter')
				tab => skyrme_eos(skyrme_matter)
				case ('neutron')
				tab => skyrme_eos(skyrme_neutron)
				case ('sym')
				tab => skyrme_eos(skyrme_sym)
				case default
				ierr = -5
				write(error_unit,*) 'unable to parse line in '//trim(model)//' for '//trim(pset)//' matter'
				return
			end select
			tab% gamma = gamma
			tab% a = a
			tab% b = b
			tab% c = c
			tab% d = d
			tab% e = e
			tab% is_loaded = .TRUE.
		end subroutine read_one
	end subroutine load_skyrme_table

end module skyrme
