program test_neutrino
    use math_lib
	use constants_def
    use constants_lib
	use superfluid_def
	use superfluid_lib
	use neutrino_def
	use neutrino_lib
    use nucchem_def
    use nucchem_lib

	type(core_neutrino_emissivity_channels) :: neu
	logical, dimension(num_core_nu_channels) :: which_channels
	integer :: sf, ierr, i, j
	real(dp) :: T, nn, np, kn, kp, Tc(max_number_sf_types)
    real(dp), parameter :: default_nuclear_radius = 1.13 ! fm
	! composition taken from HZ090 for Fe-chain accreted crust
	integer, dimension(6:14), parameter :: 	zz = [26,26,26,26,24,20,14,24,24]
    integer, dimension(6:14), parameter :: aa = [56,56,56,56,56,56,46,96,96]
	real(dp), dimension(6:14), parameter :: xn = [0.00,0.00,0.00,0.00,0.00,0.00,0.07,0.76,0.80]
	integer, dimension(2) :: Z,  N, chem_ids, charged_ids
	real(dp), dimension(2) :: Y, Yion
	integer :: ncharged
	type(composition_info_type) :: ionic
	real(dp) :: Xsum
	
    call math_init()
    call constants_init('../..','',ierr)
    if (ierr /= 0) then
        write (*,*) 'unable to initialize constants'
        stop
    end if
    
	nn = 0.16
	np = 0.02
	
	kn = (3.0*pi**2*nn)**one_third
	kp = (3.0*pi**2*np)**one_third
	
	print '(a)','setting gaps'
	call sf_startup(ierr)
	call sf_load_gaps('ns','gc','t72',ierr)
	call sf_get_results(kp,kn,Tc)
	print '(A,es11.4)','Tcns = ',Tc(neutron_1S0)
	print '(A,es11.4)','Tcps = ',Tc(proton_1S0)
	print '(A,es11.4)','Tcnt = ',Tc(neutron_3P2)
		
    print '(/,/,a)','core neutrino emissivity'
	print '(/,a)','setting emissivity to minimal cooling'
	which_channels = core_nu_minimal_cooling
	call do_one_core
	
	print '(/,a)','now adding direct Urca to minimal cooling'
	which_channels(icore_nu_dUrca) = .TRUE.
	call do_one_core
	
    print '(/,/,a)', 'crust neutrino emissivity'
    call nucchem_init(ierr)
    
	do i = 6,14
		N = [1,aa(i)-zz(i)]
		Z = [0,zz(i)]
		chem_ids = [(get_nuclide_index_from_ZN(Z(j),N(j)),j=1,2)]
		Y = [xn(i),(1.0-xn(i))/real(aa(i),dp) ]
		call compute_composition_moments(2,chem_ids,Y,ionic,Xsum,ncharged, charged_ids, Yion, exclude_neutrons=.TRUE.)
		ionic% Q = 4.0
		call do_one_crust(10.0_dp**i)
	end do
	
	call clear_composition(ionic)
    call nucchem_shutdown
	call sf_shutdown
    
	contains
	subroutine do_one_core()
		write (*,'(/,10(a11,tr1),/,10(11("="),tr1))') &
			&  'T','brem_np','brem_pp','brem_nn','mUrca_p','mUrca_n','dUrca','cooper_n','cooper_p','total'
		do i = 1,20
			T = exp10(7.0_dp+2.0_dp*(i-1.0_dp)/19.0_dp)
			call get_core_neutrino_emissivity(nn,np,T,Tc,which_channels,neu)
			write(*,'(10(es11.4,tr1))') T,neu% brem_np, neu% brem_pp, neu% brem_nn, neu% mUrca_p, &
				& neu% mUrca_n, neu% dUrca, neu% PBF_n, neu% PBF_p, neu% total
		end do
	end subroutine do_one_core
	
	subroutine do_one_crust(rho)
        use math_lib
		real(dp), intent(in) :: rho
		real(dp) :: lgr, lgT, T, nn, kn, Tc(max_number_sf_types)
		real(dp) :: chi
		type(crust_neutrino_emissivity_channels) :: eps
		integer :: ii
		
		lgr = log10(rho)
		chi = onethird*fourpi*(default_nuclear_radius*fm_to_cm)**3 * (rho*(1.0-ionic%Yn)/amu)
		nn = rho*ionic%Yn/amu/(1.0-chi)
		kn = (1.5*pi**2*nn)**onethird / cm_to_fm
		call sf_get_results(0.0_dp,kn,Tc)
		do ii = 1, 11
			lgT = 7.5 + real(ii-1,dp)/10.0_dp
			T = exp10(lgT)
			call get_crust_neutrino_emissivity(rho,T,ionic,chi,Tc(neutron_1S0),eps)
			print '(5f6.2,6es14.6)', &
			 	&   lgr,lgT,ionic%Z,ionic%A,ionic%Yn, &
			 	&   eps%total,eps%pair,eps%photo,eps%plasma,eps%brems,eps%pbf
		end do
	end subroutine do_one_crust
	
end program test_neutrino
