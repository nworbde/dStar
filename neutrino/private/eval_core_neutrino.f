module eval_core_neutrino
    use neutrino_def
    use constants_def
    real(dp), parameter :: threepithird = threepisquare**one_third

contains

    subroutine emissivity(nn,np,T,Tcs,include_channel,eps)
        use superfluid_def, only: max_number_sf_types, proton_1S0, neutron_1S0, neutron_3P2
        !INPUT: nn, np: densities (fm^-3^) of n, p
        !               T is temperature in K, Tcs are critcal temperatures
        real(dp), intent(in) :: nn,np,T,Tcs(max_number_sf_types)
        logical, dimension(num_core_nu_channels), intent(in) :: include_channel
        type(core_neutrino_emissivity_channels), intent(out) :: eps
        real(dp) :: tps, tns, tnt, vps, vns, vnt, istps, istns, T9
    
        tps = T/Tcs(proton_1S0)
        tns = T/Tcs(neutron_1S0)
        tnt = T/Tcs(neutron_3P2)
        T9 = T*1.0e-9_dp

        istps = 1.0_dp/sqrt(tps)
        istns = 1.0_dp/sqrt(tns)
        vps = 0.0_dp
        if(tps < 1.0_dp)  &
            & vps = sqrt(1.0_dp-tps)*(1.456_dp+istps*(-0.157_dp+istps*1.764_dp))
        vns = 0.0_dp
        if(tns < 1.0_dp) &
            & vns = sqrt(1.0_dp-tns)*(1.456_dp+istns*(-0.157_dp+istns*1.764_dp))
        vnt = 0.0_dp
        if(tnt < 1.0_dp) vnt = sqrt(1.0_dp-tnt)*(0.7893_dp+1.188_dp/tnt)
    
        call clear_epsilon      
        if (include_channel(icore_nu_brem)) call nubrems( &
            & nn,np,T9,vps,vns,vnt,eps% brem, &
            & eps% brem_np, eps% brem_pp, eps% brem_nn)
        
        if (include_channel(icore_nu_mUrca)) call mUrca( &
            & nn,np,T9,vps,vns,vnt,eps% mUrca, eps% mUrca_p, eps% mUrca_n)

        if (include_channel(icore_nu_dUrca)) call dUrca( &
            & nn,np,T9,vps,vns,vnt,eps% dUrca)

        if (include_channel(icore_nu_PBF)) call nupbf( &
            & nn,np,T9,vps,vns,vnt,eps% PBF, eps% PBF_n, eps% PBF_p)
    
        ! sum contributions
        eps% total = eps% brem + eps% mUrca + eps% dUrca + eps% PBF
    
        contains
        subroutine clear_epsilon()
            eps% total = 0.0_dp
            eps% brem = 0.0_dp; eps% brem_np = 0.0_dp
            eps% brem_pp = 0.0_dp; eps% brem_nn = 0.0_dp
            eps% mUrca = 0.0_dp; eps% mUrca_p = 0.0_dp; eps% mUrca_n = 0.0_dp
            eps% dUrca = 0.0_dp
            eps% PBF = 0.0_dp; eps% PBF_n = 0.0_dp; eps% PBF_p = 0.0_dp
        end subroutine clear_epsilon    
    end subroutine emissivity

    subroutine nubrems(nn,np,T9,vps,vns,vnt,eps,eps_np,eps_pp,eps_nn)
        ! bremsstrahlung
        real(dp), intent(in) :: nn,np,T9,vps,vns,vnt
        real(dp), intent(out) :: eps,eps_np, eps_pp, eps_nn
        real(dp), parameter :: Qnn = 7.4e19_dp, Qnp = 1.5e20_dp, Qpp = 7.4e19_dp
        real(dp) :: T98, rbpnp, rbppp, rbnns, rbnnt, rbnnn, rbnp, rda, rdb, rd, ap, bp, cp, ddp, cn, dn
    
        T98 = T9**8
        rbpnp = 1.0_dp
        rbppp = 1.0_dp
        if (vps > 0.0_dp) then
            ap = 0.9982 + sqrt(0.0018**2+(0.3815*vps)**2)
            bp = 0.3949 + sqrt(0.6051**2+(0.2666*vps)**2)
            cp = 0.1747 + sqrt(0.8253**2+(0.07933*vps)**2)
            ddp = 0.7333 + sqrt(0.2667**2+(0.1678*vps)**2)
            rbpnp = (1.0/2.732)*(ap*exp(1.306-sqrt(1.306**2+vps**2)) &
                & + 1.732*bp**7*exp(3.303-sqrt(3.303**2+4.0*vps**2)))
            rbppp = 0.5*(cp**2*exp(4.228-sqrt(4.228**2+(2*vps)**2)) & 
                & + ddp**7.5*exp(7.762-sqrt(7.762**2+(3.0*vps)**2)))
        end if
    
        rbnns = 1.0_dp
        if (vns > 0.0) then
            cn = 0.1747 + sqrt(0.8253**2+(0.07933*vns)**2)
            dn = 0.7333 + sqrt(0.2667**2+(0.1678*vns)**2)
            rbnns = 0.5*(cn**2*exp(4.228-sqrt(4.228**2+(2*vns)**2)) &
                & + dn**7.5*exp(7.762-sqrt(7.762**2+(3.0*vns)**2)))
        end if
    
        rbnnt = 1.0
        if (vnt > 0.0) then 
            cn = 0.1747 + sqrt(0.8253**2+(0.07933*vnt)**2)
            dn = 0.7333 + sqrt(0.2667**2+(0.1678*vnt)**2)
            rbnnt = 0.5*(cn**2*exp(4.228-sqrt(4.228**2+(2*vnt)**2)) &
                & + dn**7.5*exp(7.762-sqrt(7.762**2+(3.0*vnt)**2)))
        end if
    
        rda = 1.0
        if(vps > 0.0) rda = (0.2312 + sqrt((0.7688)**2+(0.1438*vps)**2))**5.5 &
            & *exp(3.427*(1.0-sqrt(1.0+(vps/3.427)**2)))
        rdb = 1.0
        if(vnt > 0.0) rdb = (0.2546 + sqrt((0.7454)**2 + (0.1284*vnt)**2))**5 &
            & *exp(2.701*(1.0-sqrt(1.0+(vnt/2.701)**2)))
        rd = min(rda,rdb)

        rbnnn = min(rbnns, rbnnt)
        rbnp = rbpnp
        if (rda > 0.0) rbnp = (rd/rda)*rbnp
    
        eps_np = Qnp*T98*rbnp
        eps_pp = Qpp*T98*rbppp
        eps_nn = Qnn*T98*rbnnn
        eps = T98*(Qnp*rbnp + Qpp*rbppp + Qnn*rbnnn)
    end subroutine nubrems

    subroutine mUrca(nn,np,T9,vps,vns,vnt,eps, eps_p, eps_n)
    ! modified Urca
        real(dp), intent(in) :: nn, np, T9, vps, vns, vnt
        real(dp), intent(out) :: eps, eps_p, eps_n
        real(dp), parameter :: Qn = 8.55e21_dp, Qp = 8.53e21_dp
        real(dp) :: T98, rmp, rmn, aa, ba, rmpa, rmpb, ab, bb

        T98 = T9**8
    
        ! neutron branch
        rmn = 1.0
        if (vps > 0.0) then
            aa = 0.1477 + sqrt((0.8523)**2+(0.1175*vps)**2)
            ba = 0.1477 + sqrt((0.8523)**2+(0.1297*vps)**2)
            rmn = 0.5*(aa**7.5 +ba**5.5)*exp(3.4370*(1.0-sqrt(1.0+(vps/3.4370)**2)))
        end if

        ! proton branch
        rmpa = (0.2414+sqrt((0.7586)**2+(0.1318*vps)**2))**7 &
                & *exp(5.339*(1.0-sqrt(1.0+(2.0*vps/5.339)**2)))
        rmpb = 1.0
        if (vnt > 0.) then
            ab = 0.1612+sqrt((0.8388)**2+(0.1117*vnt)**2)
            bb = 0.1612+sqrt(0.8388**2+(0.1274*vnt)**2)
            rmpb = 0.5*(ab**7+bb**5)*exp(2.398*(1.0-sqrt(1.0+(vnt/2.398)**2)))
        end if
        rmp = min(rmpa, rmpb)

        eps_n = Qn*T98*rmn 
        eps_p = Qp*T98*rmp
        eps = Qn*T98*rmn + Qp*T98*rmp
    end subroutine mUrca

    subroutine dUrca(nn,np,T9,vps,vns,vnt,eps)
    ! direct Urca
        real(dp), intent(in) :: nn, np, T9, vps, vns, vnt
        real(dp), intent(out) :: eps
        real(dp), parameter :: Q = 1.0e27
        real(dp) :: T96, rda, rdb, rd
        
        T96 = T9**6 
        rda = 1.0_dp
        if(vps > 0.0_dp)  &
            & rda = (0.2312 + sqrt((0.7688)**2 + (0.1438*vps)**2))**5.5 &
            & *exp(3.427*(1.0-sqrt(1.0+(vps/3.427)**2)))
        rdb = 1.0_dp
        if(vnt > 0.0_dp) &
            & rdb = (0.2546 + sqrt((0.7454)**2 + (0.1284*vnt)**2))**5 &
            & *exp(2.701*(1.0-sqrt(1.0+(vnt/2.701)**2)))
        rd = min(rda,rdb)
        
        eps = Q*T96*rd
    end subroutine dUrca
    
    subroutine nupbf(nn,np,T9,vps,vns,vnt,eps,eps_n,eps_p)
        real(dp), intent(in) :: nn,np,T9,vps,vns,vnt
        real(dp), intent(out) :: eps, eps_n, eps_p
        real(dp), parameter :: Q = 1.170e21_dp*3.0, a_cns=1.0, a_cnt=3.18
        real(dp) :: T97
        real(dp) :: n, x, xFn, xFp, mn, mp, fa_n, fa_ns, fa_nt, fa_p, a_cp, pFn2, PFp2, vn2, vp2, v2
        
        n = nn + np
        x = np/n
        T97 = T9**7
        
        call effective_masses(n,x,mn,mp)
        mn = mn/Mn_n
        mp = mp/Mp_n

        ! proton factor
        fa_p = 0.0_dp
        if (vps > 0.0_dp) then
            v2 = vps**2
            fa_p = (0.602*v2 + 0.5942*v2**2 + 0.288*v2**3)*sqrt(0.5547+sqrt(0.4453**2+0.0113*v2)) &
                & *exp(-sqrt(4.0*v2+2.245**2)+2.245)
        end if

        ! neutron singlet
        fa_ns = 0.0_dp
        if (vns > 0.0_dp) then
            v2 = vns**2
            fa_ns = (0.602*v2 + 0.5942*v2**2 + 0.288*v2**3)*sqrt(0.5547+sqrt(0.4453**2+0.0113*v2)) &
                & *exp(-sqrt(4.0*v2+2.245**2)+2.245)
        end if

        ! neutron triplet (only case B for now)
        fa_nt = 0.0_dp
        if (vnt > 0.0_dp) then
            v2 = vnt**2
            fa_nt = (1.204*v2+3.733*v2**2+0.3191*v2**3)/(1.0+0.3511*v2)*(0.7591+sqrt(0.2409**2+0.3145*v2))**2* &
                & exp(-sqrt(4.0*v2+0.4616**2)+0.4616)
        end if
        
        fa_n = fa_nt*a_cnt + fa_ns*a_cns

        ! v/c corrections
        xFn = threepithird*nn**one_third * hbarc_n/Mn_n
        xFp = threepithird*np**one_third * hbarc_n/Mp_n
        
        pFn2 = xFn**2
        pFp2 = xFp**2
        vn2 = pFn2/(1.0+pFn2)
        vp2 = pFp2/(1.0+pFp2)

        a_cp  = 0.0064 + 1.59*pFp2 * (mp**2 + 11.0/42.0)
        
        eps_n = Q*T97*(mn*xFn*3.0*fa_n*vn2 + mn*xFn*3.0*fa_nt)
        eps_p = Q*T97*(mp*xFp*3.0*fa_p*a_cp*vp2)
        
        eps = Q*T97*(mn*xFn*3.0*fa_n*vn2 + mn*xFn*3.0*fa_nt + mp*xFp*3.0*fa_p*a_cp*vp2)     

    contains
    subroutine effective_masses(n,x, mn,mp)
        real(dp), parameter :: fm_to_cm = 1.0e-13_dp
        real(dp), parameter :: Mneutron = 1.674927351e-24_dp
        real(dp), parameter :: Mproton  = 1.672621777e-24_dp
        
        real(dp) :: hbarc_n, Mn_n, Mp_n
        real(dp), intent(in) :: n,x
        real(dp), intent(out) :: mn, mp
        ! following constants are from APR EOS
        real(dp), parameter ::  p3 = 89.8_dp   ! MeV fm^5^
        real(dp), parameter ::  p4 = 0.457_dp  ! fm^3^
        real(dp), parameter ::  p5 = -59.0_dp  ! MeV fm^5^ 
        real(dp) :: t1, t2, t3
        
        hbarc_n = hbar*clight/mev_to_ergs/fm_to_cm
        Mn_n = Mneutron*clight**2/mev_to_ergs
        Mp_n = Mproton*clight**2/mev_to_ergs
        
        t1 = n*exp(-p4*n)
        t2 = p3+p5*(1.0-x)
        t3 = p3+x*p5
        mn = 1.0/(1.0/Mn_n + 2.0*t1*t2/hbarc_n**2)
        mp = 1.0/(1.0/Mp_n + 2.0*t1*t3/hbarc_n**2)      
    end subroutine effective_masses
    end subroutine nupbf

end module eval_core_neutrino
