program test_neutrino
	use constants_def
    use constants_lib
	use superfluid_def
	use superfluid_lib
	use core_neutrino_def
	use core_neutrino_lib

	type(core_neutrino_emissivity_channels) :: neu
	logical, dimension(num_core_nu_channels) :: which_channels
	integer :: sf, ierr, i
	real(dp) :: T, nn, np, kn, kp, Tc(max_number_sf_types)
	
    call constants_init('',ierr)
    if (ierr /= 0) then
        write (*,*) 'unable to initialize constants'
        stop
    end if
    
	nn = 0.16
	np = 0.02
	
	kn = (3.0*pi**2*nn)**one_third
	kp = (3.0*pi**2*np)**one_third
	
	print '(a)','setting gaps'
	call sf_startup('../../data',ierr)
	call sf_load_gaps('ns','gc','t72',ierr)
	call sf_get_results(kp,kn,Tc)
	print '(A,es11.4)','Tcns = ',Tc(neutron_1S0)
	print '(A,es11.4)','Tcps = ',Tc(proton_1S0)
	print '(A,es11.4)','Tcnt = ',Tc(neutron_3P2)
		
	print '(/,a)','setting emissivity to minimal cooling'
	which_channels = core_nu_minimal_cooling
	call do_one
	
	print '(/,a)','now adding direct Urca to minimal cooling'
	which_channels(icore_nu_dUrca) = .TRUE.
	call do_one
	
	call sf_shutdown
	
	contains
	subroutine do_one()
		write (*,'(/,10(a11,tr1),/,10(11("="),tr1))') &
			&  'T','brem_np','brem_pp','brem_nn','mUrca_p','mUrca_n','dUrca','cooper_n','cooper_p','total'
		do i = 1,20
			T = 10.0**(7.0+2.0*(i-1.0)/19.0)
			call get_core_neutrino_emissivity(nn,np,T,Tc,which_channels,neu)
			write(*,'(10(es11.4,tr1))') T,neu% brem_np, neu% brem_pp, neu% brem_nn, neu% mUrca_p, &
				& neu% mUrca_n, neu% dUrca, neu% PBF_n, neu% PBF_p, neu% total
		end do
	end subroutine do_one
	
	
end program test_neutrino
