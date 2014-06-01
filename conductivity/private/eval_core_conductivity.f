module eval_core_conductivity
    use conductivity_def

contains   
    function neutron_conductivity(nn,np,mneff,mpeff,T,Tcs) result (Kn)
        use constants_def
        use superfluid_def, only: max_number_sf_types, proton_1S0, neutron_1S0, neutron_3P2
        ! input nn, np: densities (fm**-3)
        ! T is temperature in K, Tcs are critical temperatures in K
        ! mneff = effective mass of neutron/neutron rest mass
        ! mpeff = effective mass of proton/proton rest mass
        real(dp), intent(in) :: nn,np,mneff,mpeff,T,Tcs(max_number_sf_types)
        real(dp) :: Kn
        real(dp) :: T8,tps,tns,tnt,istps,istns,vps,vns,vnt,yn,yp
        real(dp), parameter :: n0 = 0.16  ! fm**-3
        real(dp) :: y, Sn1, Sn2, Sp1, Sp2
        real(dp) :: kFn, kFp, kFn2, kFp2
        real(dp) :: Kn1, Kn2, Kp1, Kp2, un1, un2, up1, up2
        real(dp) :: nu_nn, nu_np
        real(dp) :: mneffinv2, mpeffinv2

        mneffinv2 = 1.0/mneff**2
        mpeffinv2 = 1.0/mpeff**2

        tps = T/Tcs(proton_1S0)
        tns = T/Tcs(neutron_1S0)
        tnt = T/Tcs(neutron_3P2)

        istps = 1./sqrt(tps)
        istns = 1./sqrt(tns)
        vps = 0.0
        if(tps < 1.0) vps = sqrt(1.-tps)*(1.456+istps*(-0.157+istps*1.764))
        vns = 0.0
        if(tns < 1.0) vns = sqrt(1.-tns)*(1.456+istns*(-0.157+istns*1.764))
        vnt = 0.0
        if(tnt < 1.0) vnt = sqrt(1.-tnt)*(0.7893+1.188/tnt)

        yn = max(vns,vnt)
        yp = vps

        kFn = (threepisquare*nn)**onethird
        kFp = (threepisquare*np)**onethird
        kFn2 = kFn**2
        kFp2 = kFp**2

        Sn1 = 14.57/kFn**1.5 * (1.0-0.0788*kFn+0.0883*kFn2)/(1.0-0.1114*kFn)
        Sn2 = 7.880/kFn2 * (1.0-0.2241*kFn+0.2006*kFn2)/(1.0-0.1742*kFn)
        Sp1 = 0.8007*kFp/kFn**2 * (1.0+31.28*kFp-0.0004285*kFp**2 &
        & + 26.85*kFn + 0.08012*kFn**2)/(1.0-0.5898*kFn+0.2368*kFn**2+0.5838*kFp**2+0.884*kFn*kFp)
        Sp2 = 0.3830*kFp**4/kFn**5.5*(1.0+102.0*kFp+53.91*kFn) &
        & /(1.0-0.7087*kFn+0.2537*kFn**2+9.404*kFp**2-1.589*kFn*kFp)

        ! in-medium corrections
        !       un1 = kFn-1.665
        ! Kn1 = mneffinv2 * (0.4583+0.892*un1**2-0.5497*un1**3-0.06205*kFp+0.04022*kFp**2+0.2122*un1*kFp)
        ! un2 = kFn-1.556
        ! Kn2 = mneffinv2 * (0.4891+1.111*un2**2-0.2283*un2**3+0.01589*kFp-0.02099*kFp**2+0.2773*un2*kFp)
        !       up1 = kFn-2.126
        ! Kp1 = mpeffinv2 * (0.04377+1.1*up1**2+0.1180*up1**3+0.1626*kFp+0.3871*up1*kFp-0.2990*up1**4)
        ! up2 = kFn-2.116
        ! Kp2 = mpeffinv2 * (0.0001313+1.248*up2**2+0.2403*up2**3 &
        !    & +0.3257*kFp+0.5536*up2*kFp-0.3237*up2**4+0.09786*up2**2*kFp)

        Kn1 = 1.0
        Kn2 = 1.0
        Kp1 = 1.0
        Kp2 = 1.0

        ! put the whole thing together
        T8 = T*1.0e-8
        nu_nn = 3.48e15*T8**2 *mneff**3 * (Sn2*Kn2*Rn2(yn)+3.0*Sn1*Kn1*(Rn1(yn)-Rn2(yn)))
        nu_np = 3.48e15*T8**2*mneff*mpeff**2* &
        & (Sp2*Kp2*Rp2(yn,yp)+0.5*Kp1*Sp1*(3.0*Rp1(yn,yp)-Rp2(yn,yp)))

        Kn = 7.2e23*T8*RC(yn)**2 /mneff * (1.0e15/(nu_nn+nu_np)) * nn/n0

    contains
        function RC(y)
            real(dp), intent(in) :: y
            real(dp) :: RC,y2
            y2 = y**2
            RC = (0.647+sqrt(0.353**2 + 0.109*y2))**1.5 * exp(1.39-sqrt(1.39**2+y2))
        end function RC
        function Rn1(y)
            real(dp), intent(in) :: y
            real(dp) :: Rn1, y2
            y2 = y**2
            Rn1 = twothird*(0.9468+sqrt(0.0532**2+0.5346*y2))**3 &
            & * exp(0.377-sqrt(0.377**2 + 4.0*y2)) &
            & + onethird*(1.0+1.351*y2)**2 &
            & * exp(0.169-sqrt(0.169**2+9.0*y2))
        end function Rn1
        function Rn2(y)
            real(dp), intent(in) :: y
            real(dp) :: Rn2,y2
            y2 = y**2
            Rn2 = 0.5*(0.6242+sqrt(0.3758**2 + 0.07198*y2))**2  &
            & * exp(3.6724 - sqrt(3.6724**2+4.0*y2))  &
            & + 0.5*(1.0+0.01211*y2)**9 &
            & * exp(7.5351-sqrt(7.5351**2+9.0*y2))
        end function Rn2
    end function neutron_conductivity

    function Rp1(yn,yp)
        real(dp), intent(in) :: yn,yp
        real(dp) :: Rp1
        real(dp) :: yn2,yp2,uplus, uminus, un, up, yminus, yplus
        yn2 = yn**2
        yp2 = yp**2
        if (yn2 > 0.0 .and. yp2 == 0.0) then
            Rp1 = (0.4459 + sqrt(0.5541**2 + 0.03016*yn2))**2 &
            & * exp(2.1178-sqrt(2.1178**2 + yn2))
        else if (yn2 == 0.0 .and. yp2 > 0.0) then
            Rp1 = 0.5*(0.3695+sqrt(0.6305**2+0.01064*yp2)) &
            & * exp(2.4451-sqrt(2.4451**2+yp2)) &
            & + 0.5*(1.0+0.1917*yp2)**1.4 &
            & * exp(4.6627-sqrt(4.6627**2+4.0*yp2))
        else if (yn2 > 0.0 .and. yp2 > 0.0) then
            yminus = min(yn, yp)
            yplus = max(yn, yp)
            uplus = sqrt(yplus**2+1.485**2)-1.485
            uminus = sqrt(yminus**2+1.485**2)-1.485
            un = sqrt(yn2+1.485**2)-1.485
            up = sqrt(yp2+1.485**2)-1.485
            Rp1 = exp(-uplus-uminus)*(0.7751+0.4823*un+0.1124*up &
            & + 0.04991*un**2+0.08513*un*up+0.01284*un**2*up) &
            & + exp(-2.0*uplus) &
            & *(0.2249+0.3539*uplus-0.2189*uminus-0.6069*un*uminus + 0.7362*up*uplus)
        else
            Rp1 = 1.0
        end if
    end function Rp1
    function Rp2(yn,yp)
        real(dp), intent(in) :: yn,yp
        real(dp) :: Rp2, yn2, yp2
        real(dp) :: uplus,uminus,un,up,yminus,yplus
        yn2 = yn**2
        yp2 = yp**2
        if (yn2>0.0 .and. yp2 == 0.0) then
            Rp2 = (0.801+sqrt(0.199**2+0.04645*yn2))**2  &
            & * exp(2.3569-sqrt(2.3569**2+yn2))
        else if (yn2 == 0.0 .and. yp2 > 0.0) then
            Rp2 = 0.0436*(sqrt(4.345**2+19.55*yp2)-3.345) &
            & * exp(2.0247-sqrt(2.0247**2 + yp2)) &
            & + 0.0654*exp(8.992-sqrt(8.992**2+1.5*yp2)) &
            & + 0.891*exp(9.627-sqrt(9.627**2+9.0*yp2))
        else if (yn2 > 0.0 .and. yp2 > 0.0) then
            yplus = max(yn,yp)
            yminus = min(yn,yp)
            uplus = sqrt(yplus**2+1.761**2)-1.761
            uminus = sqrt(yminus**2+1.761**2)-1.761
            un  = sqrt(yn2+1.761**2)-1.761
            up = sqrt(yp2+1.761**2)-1.761
            Rp2 = exp(-uplus-uminus)*(1.1032+0.8645*un+0.2042*up &
            & + 0.07937*un**2+0.1451*un*up+0.01333*un**2*up) &
            & + exp(-2.0*uplus)*(-0.1032-0.2340*uplus+0.06152*un*uplus  &
            & + 0.7533*un*uminus-1.007*up*uplus)
        else
            Rp2 = 1.0
        end if
    end function Rp2
end module eval_core_conductivity
