program test_core_eos
    use constants_def, only: dp
    use constants_lib
	use nucchem_def
	use nucchem_lib
	use superfluid_lib
	use dStar_eos_lib

	integer :: eos_handle,sf_handle,ierr
	real(dp), dimension(num_dStar_eos_results) :: res
	integer :: i
	type(core_eos_component), dimension(num_core_eos_components) :: eos_components
    logical, parameter :: dbg = .TRUE.
	real(dp) :: x
    
    call constants_init('',ierr)
	call nucchem_init('../../data',ierr)
	call sf_startup('../../data',ierr)
	call sf_load_gaps('ns','gc','t72',ierr)
	
	call dStar_eos_startup('../../data')
	eos_handle = alloc_dStar_eos_handle(ierr)
	
	call dStar_eos_set_controls(eos_handle,use_skyrme=.TRUE., &
	&	skyrme_parameter_set='s7d')
	
	x = 0.05_dp
	call write_headers
	call do_one(1.0e13_dp)
	call do_one(1.0e14_dp)
	call do_one(1.0e15_dp)

	call dStar_eos_shutdown
	call sf_shutdown
	call nucchem_shutdown
	
	contains
	subroutine write_headers()
		write (*,'(2a8,a12,10a8,/)') 'lg(rho)','lg(T)','Cv/kb/Na','lg(E)','lg(P)','lg(S)', &
			& 'chi_Rho','chi_T','G1','G3-1','del_ad','mu_e','mu_n'
	end subroutine write_headers
	
	subroutine do_one(rho)
		use constants_def, only: ln10,avogadro,boltzmann
		real(dp), intent(in) :: rho
		real(dp) :: lgr, lgT, T, Tc(max_number_sf_types)
		real(dp) :: chi,Gamma,TpT,f,u,p,s,cv,chi_rho,chi_T
		integer :: phase
		integer :: i
		
		lgr = log10(rho)
		do i = 1, 11
			lgT = 7.5 + (i-1)/10.0
			T = 10.0**lgT
			call eval_core_eos(eos_handle,rho,T,x,res,eos_components)
			write (*,'(2f8.2,es12.4,10f8.3)') lgr,lgT, &
				& res(i_Cv)/boltzmann/avogadro, &
				& res(i_lnE)/ln10,res(i_lnP)/ln10,res(i_lnS)/ln10, &
				& res(i_chiRho),res(i_chiT),res(i_Gamma1),res(i_Gamma3),  &
				& res(i_grad_ad),res(i_mu_e),res(i_mu_n)
			
                if (dbg) then
                    write (*,'(a11,8a12)') 'component: ', &
                    	& 'pressure','energy','entropy','free energy', &
                    	& 'Cv','mu','chi_rho','chi_T'
                    write (*,'(a11,8es12.4)') 'electron: ', &
                    	& eos_components(icore_eos_ele)
                    write (*,'(a11,8es12.4)') 'nucleon: ', &
                    	& eos_components(icore_eos_nucleon)
                end if
		end do
	end subroutine do_one

end program test_core_eos
