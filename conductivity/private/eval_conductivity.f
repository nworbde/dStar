module eval_conductivity
    use conductivity_def

contains
    
    subroutine conductivity(rho,T,chi,Gamma,eta,mu_e,ionic,Tns, &
    &   kappa,which_ee,which_eQ,K_components)
        use iso_fortran_env, only : error_unit
        use constants_def
        use nucchem_def, only: composition_info_type
        use neutron_conductivity
        use electron_conductivity
        use photon_conductivity
        
        real(dp), intent(in) :: rho, T, chi, Gamma, eta, mu_e   ! mu_e is in MeV
        type(composition_info_type), intent(in) :: ionic
        real(dp), intent(in) :: Tns ! superfluid critical temperature, neutrons
        type(conductivity_components), intent(out) :: kappa
        integer, intent(in) :: which_ee, which_eQ
        logical, dimension(num_conductivity_channels), intent(in) :: &
        &   K_components
        real(dp) :: nn, nion, ne, kF, xF, eF, Gamma_e, kappa_e_pre
        real(dp) :: kappa_n_pre, kfn, mnstar, tau, v, R
        real(dp) :: nu, nu_c, K_opacity
        integer :: ierr

        call clear_kappa
        
        ! electrons
        if (any(K_components(icond_ee:icond_eQ))) then
            ne = rho/amu*ionic%Ye
            kF = (threepisquare*ne)**onethird
            xF = hbar*kF/Melectron/clight
            eF = sqrt(1.0_dp+xF**2)
            Gamma_e = Gamma/ionic%Z53
            kappa_e_pre = onethird*pi**2*boltzmann**2*T*ne/Melectron/eF
            nu = 0.0_dp
            nu_c = 0.0_dp
            if (K_components(icond_ee)) then
                select case (which_ee)
                case(icond_sy06)
                    nu_c = ee_SY06(ne,T)
                case(icond_pcy)
                    nu_c = ee_PCY(ne,T)
                case default
                    write(error_unit,*)  &
                    &   'unknown option for ee scattering: using SY06'
                    nu_c = ee_SY06(ne,T)
                end select
                nu = nu + nu_c
                kappa% ee = kappa_e_pre/nu_c
            end if
            if (K_components(icond_ei)) then
                nu_c = eion( &
                &   kF,Gamma_e,eta,ionic%Ye,ionic%Z,ionic%Z2,ionic%Z53,ionic%A)
                nu = nu + nu_c
                if (nu_c > tiny(1.0_dp)) kappa% ei = kappa_e_pre/nu_c
            end if
            if (K_components(icond_eQ)) then
                select case(which_eQ)
                case(icond_eQ_potekhin)
                    nu_c = eQ(kF,T,ionic%Ye,ionic%Z,ionic%Z2,ionic%A,ionic%Q)
                case(icond_eQ_page)
                    nu_c = eQ_page( &
                    &    kF,T,ionic%Ye,ionic%Z,ionic%Z2,ionic%A,ionic%Q)
                end select
                nu = nu + nu_c
                if (ionic%Q > 1.0e-8_dp) kappa% eQ = kappa_e_pre/nu_c
            end if
            kappa% electron_total = kappa_e_pre/nu
        end if

        ! neutrons
        if (any(K_components(icond_sf:icond_nQ)) .and. ionic% Yn > 0) then
            nu = 0.0_dp
            nu_c = 0.0_dp
            nn = rho/amu*ionic% Yn/(1.0-chi)
            kFn = (threepisquare*nn)**onethird
            nion = rho/amu*(1.0-ionic%Yn)/ionic%A   
            mnstar = Mneutron
            kappa_n_pre = onethird*pi**2*boltzmann**2*T*nn/mnstar
            
            ! superfluid reduction factor
            R = 1.0
            if (T < Tns) then
                tau = T/Tns
                v = sqrt(1.0-tau)*(1.456-0.157/sqrt(tau) + 1.764/tau)
                R = (0.4186+sqrt(1.007**2+(0.5010*v)**2))**2.5 * &
                &   exp(1.456-sqrt(1.456**2+v**2))
                if (R < 1.0E-10) R = 0.0
            end if
            
            if (K_components(icond_nQ)) then
                nu_c = n_imp(nn,nion,T,ionic,ierr)
                kappa% nQ = kappa_n_pre/nu_c*R
                nu = nu + nu_c
            end if
            if (K_components(icond_nn)) then
                nu_c = n_phonon(nn,nion,T,ionic, ierr)
                kappa% np = kappa_n_pre/nu_c*R
                nu = nu + nu_c
            end if
            kappa% neutron_total = kappa_n_pre/nu*R

            ! superfluid phonons
            if (K_components(icond_sf) .and. T < Tns) then
                kappa% sf =  sPh(nn/density_n,nion/density_n,T,ionic)
            end if
        end if

        ! photons
        K_opacity = 0.0_dp
        if (K_components(icond_kap)) then
            kappa% kap = Rosseland_kappa(rho,T,mu_e,ionic)
            K_opacity = 4.0_dp*onethird*arad*clight*T**3/rho/kappa% kap
        end if
        
        kappa% total = kappa% electron_total  &
        &   + kappa% neutron_total + kappa% sf &
        &   + K_opacity

        contains
        subroutine clear_kappa()
            kappa% total = 0.0_dp
            kappa% ee  = 0.0_dp
            kappa% ei = 0.0_dp
            kappa% eQ = 0.0_dp
            kappa% sf = 0.0_dp
            kappa% kap = 0.0_dp
            kappa% nQ = 0.0_dp
            kappa% np = 0.0_dp
            kappa% electron_total = 0.0_dp
            kappa% neutron_total = 0.0_dp
        end subroutine clear_kappa
    end subroutine conductivity

end module eval_conductivity
