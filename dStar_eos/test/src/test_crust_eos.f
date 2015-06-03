program test_crust_eos
    use constants_def, only: dp
    use constants_lib
	use nucchem_def
	use nucchem_lib
	use superfluid_lib
	use dStar_eos_lib

	integer :: eos_handle,sf_handle,ierr
	integer, dimension(3) :: Z,A,N,chem_ids, charged_ids
	real(dp), dimension(3) :: Y, Yion
	real(dp), dimension(num_dStar_eos_results) :: res
	type(composition_info_type) :: ionic
	integer :: i, ncharged
	real(dp) :: Xsum
	type(crust_eos_component), dimension(num_crust_eos_components) :: eos_components
    logical, parameter :: dbg = .FALSE.
    
	! n/o16/fe56 composition
	Z = [0.0,8.0,26.0]; A = [1.0,16.0,56.0]
	N = A-Z
	
    call constants_init('',ierr)
	call nucchem_init('../../data',ierr)
	call sf_startup('../../data',ierr)
	call sf_load_gaps('ns','gc','t72',ierr)
	
	call dStar_eos_startup('../../data')
	eos_handle = alloc_dStar_eos_handle(ierr)
		
	! get indices
	chem_ids = [(get_nuclide_index_from_ZN(Z(i),N(i)),i=1,3)]

	! O-Fe mixture
	Y = [0.0,0.5,0.5]/real(A,dp)
	call compute_composition_moments(3,chem_ids,Y,ionic,Xsum,ncharged, charged_ids, Yion,  &
			& exclude_neutrons = .TRUE.)
	call write_headers
	call do_one(1.0e2_dp)
	call do_one(1.0e5_dp)
	call do_one(1.0e7_dp)
	call do_one(1.0e9_dp)
	call do_one(1.0e11_dp)

	! add neutrons
	Y = [0.5,0.0,0.5]/real(A,dp)
	call compute_composition_moments(3,chem_ids,Y,ionic,Xsum,ncharged, charged_ids, Yion,  &
		& exclude_neutrons=.TRUE.)	
	call write_headers
	call do_one(1.0e11_dp)
	call do_one(1.0e13_dp)
	
	! extreme neutron-rich composition
	Y = [0.95,0.0,0.05]/real(A,dp)
	call compute_composition_moments(3,chem_ids,Y,ionic,Xsum,ncharged, charged_ids, Yion,  &
		& exclude_neutrons=.TRUE.)	
	call write_headers
	call do_one(1.0e12_dp)
	call do_one(1.0e14_dp)

    ! make trouble
	! n/o16/fe56 composition
!     Z = [0.0,12.0,26.0]; A = [1.0,48.0,56.0]
!     N = A-Z
!     chem_ids = [(get_nuclide_index_from_ZN(Z(i),N(i)),i=1,3)]
!     Y = [0.8,0.2,0.0]/real(A,dp)
!     call compute_composition_moments(3,chem_ids,Y,ionic,Xsum,ncharged, charged_ids, Yion,  &
!         & exclude_neutrons=.TRUE.)
!     call write_headers
!     call do_one(1.0e12_dp)
!     call do_one(1.0e14_dp)

	call sf_shutdown
	call nucchem_shutdown
	
	contains
	subroutine write_headers()
		write (*,'(/,"<Z> =",f7.2,"; <A> =",f7.2)') ionic% Z, ionic% A
		write (*,'(2a8,3a12,10a8,/)') 'lg(rho)','lg(T)','Gamma','Theta','Cv/kb/Na','lg(E)','lg(P)','lg(S)', &
			& 'chi_Rho','chi_T','G1','G3-1','del_ad','mu_e','mu_n'
	end subroutine write_headers
	
	subroutine do_one(rho)
		use constants_def, only: ln10,avogadro,boltzmann
		real(dp), intent(in) :: rho
		real(dp) :: lgr, lgT, T, Tc(max_number_sf_types)
		real(dp) :: chi,Gamma,TpT,f,u,p,s,cv,chi_rho,chi_T,chk
		integer :: phase
		integer :: i
		
		chi = use_default_nuclear_size
		lgr = log10(rho)
		do i = 1, 11
			lgT = 7.5 + (i-1)/10.0
			T = 10.0**lgT
			call eval_crust_eos(eos_handle,rho,T,ionic,ncharged,charged_ids,Yion, &
				& res, phase, chi, eos_components)
			chk = exp(res(i_lnP)-res(i_lnE))/rho
			write (*,'(2f8.2,3es12.4,10f8.3)')  &
				& lgr,lgT,res(i_Gamma),res(i_Theta),res(i_Cv)/boltzmann/avogadro,res(i_lnE)/ln10,res(i_lnP)/ln10,res(i_lnS)/ln10, &
				& res(i_chiRho),res(i_chiT),res(i_Gamma1),res(i_Gamma3),  &
				& res(i_grad_ad), res(i_mu_e),res(i_mu_n)
			
                if (dbg) then
                    write (*,'(a11,7a12)') 'component: ','pressure','energy','entropy','free energy','Cv','chi_rho','chi_T'
                    write (*,'(a11,7es12.4)') 'electron: ',eos_components(icrust_eos_ep)
                    write (*,'(a11,7es12.4)') 'ion: ',eos_components(icrust_eos_ion)
                    write (*,'(a11,7es12.4)') 'neutron: ',eos_components(icrust_eos_neutron)
                    write (*,'(a11,7es12.4)') 'radiation: ',eos_components(icrust_eos_radiation)
                end if
		end do
	end subroutine do_one

end program test_crust_eos
