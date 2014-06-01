module crust_eos_lib
	use crust_eos_def
	
	contains
	subroutine crust_eos_startup(data_dir)
		use helm_alloc
		character(len=*), intent(in) :: data_dir
		integer, parameter :: imax = 261, jmax = 101  ! dimensions of our version of helm table
		integer :: ierr
		call crust_eos_def_init
		call alloc_helm_table(eos_ht, imax, jmax, ierr)
		if (ierr /= 0) then
			write (*,*) 'unable to alloc helm table'
			return
		end if
		call read_helm_table(eos_ht,data_dir,ierr)
		if (ierr /=0) then
			write (*,*) 'unable to read helm table'
			return
		end if
	end subroutine crust_eos_startup

	function alloc_crust_eos_handle(ierr)
		integer, intent(out) :: ierr
		integer :: alloc_crust_eos_handle
		alloc_crust_eos_handle = do_alloc_crust_eos(ierr)
	end function alloc_crust_eos_handle

	subroutine free_crust_eos_handle(handle)
		integer, intent(in) :: handle
		call do_free_crust_eos(handle)
	end subroutine free_crust_eos_handle

	subroutine crust_eos_ptr(handle, rq, ierr)
		integer, intent(in) :: handle
		type(crust_eos_general_info), pointer :: rq
		integer, intent(out) :: ierr
		call get_crust_eos_ptr(handle,rq,ierr)
	end subroutine crust_eos_ptr
	
	subroutine crust_eos_set_gamma_melt(handle,melt_pt,ierr)
		use alert_lib
		integer, intent(in) :: handle
		real, intent(in) :: melt_pt
		integer, intent(out) :: ierr
		type(crust_eos_general_info), pointer :: rq
		
		call crust_eos_ptr(handle,rq,ierr)
		if (ierr /= 0) return
		if (.not.rq% in_use) then
			ierr = -3
			call alert(ierr, 'bad handle: did you forget to initialize?')
			return
		end if
		rq% Gamma_melt = melt_pt
	end subroutine crust_eos_set_gamma_melt

	subroutine crust_eos_set_rsi_melt(handle,melt_pt,ierr)
		use alert_lib
		integer, intent(in) :: handle
		real, intent(in) :: melt_pt
		integer, intent(out) :: ierr
		type(crust_eos_general_info), pointer :: rq
	
		call crust_eos_ptr(handle,rq,ierr)
		if (ierr /= 0) return
		if (.not.rq% in_use) then
			ierr = -3
			call alert(ierr, 'bad handle: did you forget to initialize?')
			return
		end if
		rq% rsi_melt = melt_pt
	end subroutine crust_eos_set_rsi_melt

	subroutine crust_eos_set_abun_threshold(handle,Ythresh,ierr)
		use alert_lib
		integer, intent(in) :: handle
		real, intent(in) :: Ythresh
		integer, intent(out) :: ierr
		type(crust_eos_general_info), pointer :: rq

		call crust_eos_ptr(handle,rq,ierr)
		if (ierr /= 0) return
		if (.not.rq% in_use) then
			ierr = -3
			call alert(ierr, 'bad handle: did you forget to initialize?')
			return
		end if
		rq% Ythresh = Ythresh
	end subroutine crust_eos_set_abun_threshold
	
	subroutine eval_crust_eos( &
			& crust_eos_handle,rho,T,ionic,ncharged,charged_ids,Yion, &
			& res,phase,chi,components)
		use alert_lib
		use nucchem_def, only: composition_info_type
		use electron_eos
		use ion_eos
		use neutron_eos
		use phys_constants
		
		integer, intent(in) :: crust_eos_handle
		real, intent(in) :: rho,T
		type(composition_info_type), intent(in) :: ionic
		integer, intent(in) :: ncharged
		integer, dimension(ncharged), intent(in) :: charged_ids
			! ids of the charged species
		real, dimension(ncharged), intent(in) :: Yion
			! renormalized abunances of charged species Yion = Y/(1-Yn)
		real, dimension(num_crust_eos_results) :: res
		integer, intent(out) :: phase
		real, intent(inout) :: chi
		type(crust_eos_component), intent(out), dimension(num_crust_eos_components), optional :: components
			! volume fraction of nucleus; if input with value use_default_nuclear_size, chi is computed
			! and the new value is returned in res.  Otherwise, the value that is input is used by the code and 
			! unaltered.
		type(crust_eos_general_info), pointer :: rq		
		real :: ne,rs,Gamma_e,nek, nekT,f_e, u_e, p_e, s_e, cv_e, dpr_e, dpt_e, eta_e,mu_e
		real :: f_ex, u_ex, p_ex, s_ex, cv_ex, dpr_ex, dpt_ex, uexfac, pexfac, sexfac
		real :: n_i,nik,nikT,f_i,u_i,p_i,s_i,cv_i,dpr_i,dpt_i,uifac,pifac,sifac
		real :: nn,f_n,u_n,p_n,s_n,cv_n,dpr_n,dpt_n,mu_n,unfac,pnfac,snfac
		real :: Gamma,ionQ,p,u,s,cv,dpr,dpt,gamma3m1,gamma1,grad_ad,cp
		integer :: ierr
		
		call crust_eos_ptr(crust_eos_handle, rq, ierr)
		if (ierr /= 0) return
		
		! electrons... some ion quantities are defined in terms of these as well.
		ne = rho*ionic%Ye/amu
		rs = (3.0/fourpi/ne)**onethird / a_Bohr
		Gamma_e = 2.0*Rydberg/boltzmann/T/rs
		
		Gamma = 0.0
		ionQ = 0.0
				
		! set the nuclear size
		if (chi == use_default_nuclear_size) then
			chi = onethird*fourpi*(default_nuclear_radius*fm_to_cm)**3 * (rho*(1.0-ionic%Yn)/amu)
		end if
		
		! electrons
		call get_helm_eos_results(rho,T,ionic,f_e,u_e,p_e,s_e,cv_e,dpr_e,dpt_e,eta_e)
		call ee_exchange(Gamma_e, rs, f_ex, u_ex, p_ex, s_ex, cv_ex, dpr_ex, dpt_ex)
		nek = ne*boltzmann
		nekT = nek*T
		uexfac = nekT/rho
		pexfac = nekT
		sexfac = nek/rho
		! u = u + u_e + u_ex*nekT/rho
		! p = p + p_e + p_ex*nekT
		! s = s + s_e + s_ex*nek/rho
		! cv = cv + cv_e + cv_ex*nek/rho
		! dpr = dpr +dpr_e + dpr_ex*nekT
		! dpt = dpt +dpt_e + dpt_ex*nekT
!		u_ex = 0.0; p_ex = 0.0; s_ex = 0.0; cv_ex = 0.0; dpr_ex = 0.0; dpt_ex = 0.0
		
		! ions
		call ion_mixture(rq,rs,Gamma_e,ionic,ncharged,charged_ids,Yion, &
			& Gamma,ionQ,phase,f_i,u_i,p_i,s_i,cv_i,dpr_i,dpt_i,ierr)
		if (ierr /= 0 .and. ierr /= strong_quantum_effects) write(*,'(a)') alert_message
		n_i = ne/ionic%Z
		nik = n_i*boltzmann
		nikT = nik*T
		uifac = nikT/rho
		pifac = nikT
		sifac = nik/rho
		! u = u + u_i*nikT/rho
		! p = p + p_i*nikT
		! s = s + s_i*nik/rho
		! cv = cv + cv_i*nik/rho
		! dpr = dpr + dpr_i*nikT
		! dpt = dpt_i*nikT
!		p_i = 0.0; s_i = 0.0; u_i = 0.0; cv_i = 0.0; dpr_i = 0.0; dpt_i = 0.0
		
		! neutrons
		nn = rho*ionic%Yn/amu/(1.0-chi)
			! local density of neutrons
		! kn = (0.5*threepisquare*nn)**onethird / cm_to_fm
		! call sf_get_results(0.0,kn,Tc)
		call MB77(nn,T,f_n,u_n,p_n,s_n,cv_n,dpr_n,dpt_n)
		unfac = avogadro*ionic% Yn
		pnfac = 1.0
		snfac = unfac
		
!		p_n = 0.0; s_n = 0.0; u_n = 0.0; cv_n = 0.0; dpr_n = 0.0; dpt_n = 0.0
		! u = u + u_n*ionic% Yn*avogadro	! here we are getting the mean energy per unit mass for the mixture.
		! p = p + p_n
		! s = s + s_n*ionic% Yn*avogadro
		! cv = cv + cv_n*ionic% Yn*avogadro
		! dpr = dpr + dpr_n
		! mu_n = u_n - T*s_n + p_n
		
		! stuff results into output structure
		p = p_e + p_ex*pexfac + p_i*pifac + p_n*pnfac
		u = u_e + u_ex*uexfac + u_i*uifac + u_n*unfac
		s = s_e + s_ex*sexfac + s_i*sifac + s_n*snfac
		cv = cv_e + cv_ex*sexfac + cv_i*sifac + cv_n*snfac
		dpr = dpr_e + dpr_ex*pexfac + dpr_i*pifac + dpr_n*pnfac
		dpt = dpt_e + dpt_ex*pexfac + dpt_i*pifac + dpt_n*pnfac
		gamma3m1 = dpt/rho/T/cv
		gamma1 = (dpr + dpt*gamma3m1)/p
		grad_ad = gamma3m1/gamma1
		cp = gamma1*p*cv/dpr
		mu_e = (eta_e * boltzmann*T)*ergs_to_mev
		if (ionic% Yn < rq% Ythresh) then
			mu_n = 0.0
		else
			mu_n = (u_n -T*s_n + p_n/nn)*ergs_to_mev
		end if
		
		res(i_lnP) = log(p)
		res(i_lnE) = log(u)
		res(i_lnS)	= log(s)
		res(i_grad_ad) = grad_ad
		res(i_chiRho) = dpr/p
		res(i_chiT) = dpt/p
		res(i_Cp) = cp
		res(i_Cv) = cv
		res(i_Chi) = chi
		res(i_Gamma) = Gamma
		res(i_Theta) = 1.0/ionQ
		res(i_mu_e) = mu_e
		res(i_mu_n) = mu_n
		res(i_Gamma1) = gamma1
		res(i_Gamma3) = gamma3m1
		
		if (present(components)) then
			components(icrust_eos_ep)% P = p_e + p_ex*pexfac
			components(icrust_eos_ep)% E = u_e + u_ex*uexfac
			components(icrust_eos_ep)% S = s_e + s_ex*sexfac
			components(icrust_eos_ep)% F = f_e + f_ex*uexfac
			components(icrust_eos_ep)% Cv = cv_e + cv_ex*sexfac
			components(icrust_eos_ep)% dPdlnRho = (dpr_e + dpr_ex*pexfac)/(p_e+p_ex*pexfac)
			components(icrust_eos_ep)% dPdlnT = (dpt_e + dpt_ex*pexfac)/(p_e+p_ex*pexfac)

			components(icrust_eos_ion)% P = p_i*pifac
			components(icrust_eos_ion)% E = u_i*uifac
			components(icrust_eos_ion)% S = s_i*sifac
			components(icrust_eos_ion)% F = f_i*uifac
			components(icrust_eos_ion)% Cv = cv_i*sifac
			components(icrust_eos_ion)% dPdlnRho = dpr_i/p_i
			components(icrust_eos_ion)% dPdlnT = dpt_i/p_i

			if (ionic% Yn < rq% Ythresh) then
				components(icrust_eos_neutron) = crust_eos_component(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
			else
				components(icrust_eos_neutron)% P = p_n*pnfac
				components(icrust_eos_neutron)% E = u_n*unfac
				components(icrust_eos_neutron)% S = s_n*snfac
				components(icrust_eos_neutron)% F = f_n*unfac
				components(icrust_eos_neutron)% Cv = cv_n*snfac
				components(icrust_eos_neutron)% dPdlnRho = dpr_n/p_n
				components(icrust_eos_neutron)% dPdlnT = dpt_n/p_n
			end if
		end if
	end subroutine eval_crust_eos

end module crust_eos_lib
