!..here is the tabular helmholtz free energy eos:
!..
!..routine helmeos computes the pressure, energy and entropy via tables

      module helm
      implicit none
      
      contains


      subroutine helmeos2(
     >         T, logT, Rho, logRho, Zfrac, Xfrac, abar_in, zbar_in, helm_res, 
     >         clip_to_table_boundaries, ierr)
      use dStar_eos_def
      use const_def, only: pi
      implicit none
      double precision, intent(in) :: T, logT, Rho, logRho
      double precision, intent(in) :: Zfrac, Xfrac, abar_in, zbar_in
      double precision, intent(out) :: helm_res(num_helm_results)
      logical, intent(in) :: clip_to_table_boundaries
      integer, intent(out) :: ierr
      logical :: skip_elec_pos
      
      !double precision, parameter :: logRho1 = -10.0d0
      !double precision, parameter :: logRho2 = -10.5d0
      !double precision, parameter :: logT1 = 4.8d0
      !double precision, parameter :: logT2 = 4.5d0

      double precision, parameter :: logT1 = 5.0d0
      double precision, parameter :: logT2 = 4.5d0
      double precision, parameter :: logQ1 = 4d0
      double precision, parameter :: logQ2 = 3d0
      double precision, parameter :: logRho1 = logQ1 + 2*logT2 - 12
      double precision, parameter :: logRho2 = logQ2 + 2*logT2 - 12
      
      double precision :: dx, dy, dlogT2, dlogQ2, dist, alfa, beta, logQ, P, x
      double precision, dimension(num_helm_results) :: helm_res_alfa, helm_res_beta
      
      logical, parameter :: dbg = .false.
      
      include 'formats.dek'
      
      ! alfa = 0 for with ele_pos,
      ! alfa = 1 for without ele_pos,
      ! otherwise, blend.
      
      logQ = logRho - 2*logT + 12
      
      if (logT >= logT1) then ! above transition
         if (dbg) write(*,*) 'logT >= logT1'
         alfa = 0 ! full on
      else if (logT >= logT2) then ! in temperature transition region
         if (dbg) then
            write(*,*) 'logT >= logT2'
            write(*,1) 'logRho1', logRho1
            write(*,1) 'logRho2', logRho2
            write(*,1) 'logRho', logRho
            write(*,1) 'logQ', logQ
         end if
         if (logQ >= logQ1) then
            alfa = 0 ! full on
         else if (logQ < logQ2) then ! upper edge of region
            alfa = (logT - logT1) / (logT2 - logT1)
         else ! corner
            dlogT2 = ((logT - logT2) / (logT1 - logT2))**2
            dlogQ2 = (logQ - logQ2) / (logQ1 - logQ2)**2
            dist = sqrt(dlogT2 + dlogQ2) ! dist from (Q2,T2) corner
            alfa = max(0d0, 1d0 - dist)
         end if
      else ! logT < logT2
         if (dbg) write(*,*) 'logT < logT2'
         if (logRho >= logRho1) then
            alfa = 0 ! full on
         else if (logRho > logRho2) then
            alfa = (logRho - logRho1) / (logRho2 - logRho1)
         else
            alfa = 1 ! full off
         end if


      end if
      
      if (alfa > 0 .and. alfa < 1) alfa = 0.5d0*(1 - cos(pi*alfa))
      beta = 1 - alfa
      if (dbg) write(*,1) 'HELM elect-pos: alfa, beta', alfa, beta
      
      if (beta > 0) then ! eval with ele_pos
         skip_elec_pos = .false.
         call helmeos2aux(
     >         T, logT, Rho, logRho, Zfrac, Xfrac, abar_in, zbar_in, helm_res_beta, 
     >         clip_to_table_boundaries, skip_elec_pos, ierr)
         if (ierr /= 0 .or. helm_res_beta(h_stot) <= 0) then
            if (dbg) then
               write(*,1) 'T', T
               write(*,1) 'logT', logT
               write(*,1) 'Rho', Rho
               write(*,1) 'logRho', logRho
               write(*,1) 'abar', abar_in
               write(*,1) 'zbar', zbar_in
               write(*,*) 'with ele pos'
               write(*,*)
               write(*,1) 'stot', helm_res_beta(h_stot)
               write(*,1) 'sgas', helm_res_beta(h_sgas)
               write(*,1) 'srad', helm_res_beta(h_srad)
               write(*,1) 'sion', helm_res_beta(h_sion)
               write(*,1) 'sele', helm_res_beta(h_sele)
               write(*,1) 'scoul', helm_res_beta(h_scou)
               write(*,*)
               stop
            end if
            ierr = 0
            alfa = 1
            beta = 0
         end if
      end if
      
      if (alfa > 0) then ! eval without ele_pos
         skip_elec_pos = .true.
         call helmeos2aux(
     >         T, logT, Rho, logRho, Zfrac, Xfrac, abar_in, zbar_in, helm_res_alfa, 
     >         clip_to_table_boundaries, skip_elec_pos, ierr)
         if (ierr /= 0) return
         if (helm_res_alfa(h_stot) <= 0) then
            ierr = -1
            if (dbg) write(*,1) 'without ele_pos, helm_res_alfa(h_stot)', helm_res_alfa(h_stot)
            return
         end if
      end if
      
      if (alfa == 1) then
         helm_res = helm_res_alfa
         return
      end if
      
      if (beta == 1) then
         helm_res = helm_res_beta
         return
      end if

      helm_res = alfa*helm_res_alfa + beta*helm_res_beta
      ! redo the gammas, etc. to preserve consistency
      P = helm_res(h_ptot)
      if (dbg) then
         write(*,1) 'lgP blend', log10(P)
         write(*,1) 'lgP with', log10(helm_res_beta(h_ptot))
         write(*,1) 'lgP skip', log10(helm_res_alfa(h_ptot))
         write(*,*)
      end if
      helm_res(h_chit) = helm_res(h_dpt)*T/P
      helm_res(h_chid) = helm_res(h_dpd)*rho/P
      x = helm_res(h_dpt)/(helm_res(h_det)*rho)
      helm_res(h_gam3) = 1d0 + x
      helm_res(h_gam1) = helm_res(h_chit)*x + helm_res(h_chid)
      helm_res(h_nabad) = x/helm_res(h_gam1)
      helm_res(h_cp) = helm_res(h_cv)*helm_res(h_gam1)/helm_res(h_chid)
     
      end subroutine helmeos2


      subroutine helmeos2aux(
     >         temp_in, logtemp_in, den_in, logden_in, Zfrac, Xfrac, abar_in, zbar_in, helm_res, 
     >         clip_to_table_boundaries, must_skip_elec_pos, ierr)

      use dStar_eos_def
      use const_def
      use utils_lib, only: is_bad_num
      
      implicit none

      double precision, intent(in) :: temp_in, logtemp_in, den_in, logden_in
      double precision, intent(in) :: Zfrac, Xfrac, abar_in, zbar_in
      double precision, intent(out) :: helm_res(num_helm_results)
      logical, intent(in) :: clip_to_table_boundaries, must_skip_elec_pos
      integer, intent(out) :: ierr
      
      double precision :: Am, Zm, Yfrac, dabar_dlnY, dzbar_dlnY
      double precision :: dabar_dlnY_X, dzbar_dlnY_X, dabar_dlnY_Z, dzbar_dlnY_Z
      double precision :: h ! = planck_h
      type (Helm_Table), pointer :: ht

!..declare local variables
      include 'helm_declare_local_variables.dek'
      

!..given a temperature temp [K], density den [g/cm**3], and a composition 
!..characterized by abar and zbar, this routine returns most of the other 
!..thermodynamic quantities. of prime interest is the pressure [erg/cm**3], 
!..specific thermal energy [erg/gr], the entropy [erg/g/K], along with 
!..their derivatives with respect to temperature, density, abar, and zbar.
!..other quantites such the normalized chemical potential eta (plus its
!..derivatives), number density of electrons and positron pair (along 
!..with their derivatives), adiabatic indices, specific heats, and 
!..relativistically correct sound speed are also returned.
!..
!..this routine assumes planckian photons, an ideal gas of ions, 
!..and an electron-positron gas with an arbitrary degree of relativity
!..and degeneracy. interpolation in a table of the helmholtz free energy
!..is used to return the electron-positron thermodynamic quantities.
!..all other derivatives are analytic.
!..
!..references: cox & giuli chapter 24 ; timmes & swesty apj 1999

!..this routine assumes a call to subroutine read_helm_table has
!..been performed prior to calling this routine.


!..declare

      double precision abar, zbar, temp, logtemp, den, logden
      logical skip_elec_pos
      
!..for the interpolations
      integer          iat, jat
      double precision dth, dt2, dti, dt2i, dt3i, dd, dd2, ddi, dd2i, dd3i, 
     1                 xt, xd, mxt, mxd, fi(36), 
     2                 din, dindd, dinda, dindz, dindda, dinddz, dindaa, 
     3                 dindaz, dindzz, dinddaa, dinddaz, 
     2                 w0t, w1t, w2t, w0mt, w1mt, w2mt, 
     3                 w0d, w1d, w2d, w0md, w1md, w2md, 
     4                 dpepdd_in, dpepddd_in, dpepddt_in

      double precision psi0, dpsi0, ddpsi0, dddpsi0, 
     1                 psi1, dpsi1, ddpsi1, dddpsi1, 
     2                 psi2, dpsi2, ddpsi2, dddpsi2, 
     3                 h5

      double precision xpsi0, xdpsi0, xddpsi0, 
     1                 xpsi1, xdpsi1, xddpsi1, h3

      double precision si0t, si1t, si2t, si0mt, si1mt, si2mt, 
     1                 si0d, si1d, si2d, si0md, si1md, si2md, 
     2                 dsi0t, dsi1t, dsi2t, dsi0mt, dsi1mt, dsi2mt, 
     3                 dsi0d, dsi1d, dsi2d, dsi0md, dsi1md, dsi2md, 
     4                 ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, 
     5                 ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md, 
     6                 dddsi0t, dddsi1t, dddsi2t, 
     7                 dddsi0mt, dddsi1mt, dddsi2mt, 
     8                 dddsi0d, dddsi1d, dddsi2d, 
     9                 dddsi0md, dddsi1md, dddsi2md

      double precision free, df_d, df_t, df_dd, df_tt, df_dt, 
     1                 df_ttt, df_dtt, df_ddt, df_ddd



!..quintic hermite polynomial statement functions
!..psi0 and its derivatives
      psi0(z)    = z**3 * ( z * (-6.0d0*z + 15.0d0) - 10.0d0) + 1.0d0
      dpsi0(z)   = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
      ddpsi0(z)  = z* ( z*( -120.0d0*z + 180.0d0) - 60.0d0)
      dddpsi0(z) = z*( -360.0d0*z + 360.0d0) - 60.0d0


!..psi1 and its derivatives
      psi1(z)    = z* (z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
      dpsi1(z)   = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
      ddpsi1(z)  = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)
      dddpsi1(z) = z * (-180.0d0*z + 192.0d0) - 36.0d0


!..psi2  and its derivatives
      psi2(z)    = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
      dpsi2(z)   = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
      ddpsi2(z)  = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) +2.0d0)
      dddpsi2(z) = 0.5d0*(z * (-60.0d0*z + 72.0d0) - 18.0d0)


!..biquintic hermite polynomial statement function
      h5(i, j, w0t, w1t, w2t, w0mt, w1mt, w2mt, w0d, w1d, w2d, w0md, w1md, w2md)=
     1       fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t
     2     + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt
     3     + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t
     4     + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt
     5     + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t
     6     + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt
     7     + fi(13) *w1d*w0t   + fi(14) *w1md*w0t
     8     + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt
     9     + fi(17) *w2d*w0t   + fi(18) *w2md*w0t
     &     + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt
     1     + fi(21) *w1d*w1t   + fi(22) *w1md*w1t
     2     + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt
     3     + fi(25) *w2d*w1t   + fi(26) *w2md*w1t
     4     + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt
     5     + fi(29) *w1d*w2t   + fi(30) *w1md*w2t
     6     + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt
     7     + fi(33) *w2d*w2t   + fi(34) *w2md*w2t
     8     + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt



!..cubic hermite polynomial statement functions
!..psi0 & derivatives
      xpsi0(z)   = z * z * (2.0d0*z - 3.0d0) + 1.0
      xdpsi0(z)  = z * (6.0d0*z - 6.0d0)
      xddpsi0(z) = 12.0d0*z - 6.0d0


!..psi1 & derivatives
      xpsi1(z)   = z * ( z * (z - 2.0d0) + 1.0d0)
      xdpsi1(z)  = z * (3.0d0*z - 4.0d0) + 1.0d0
      xddpsi1(z) = 6.0d0*z - 4.0d0


!..bicubic hermite polynomial statement function
      h3(i, j, w0t, w1t, w0mt, w1mt, w0d, w1d, w0md, w1md) = 
     1       fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t 
     2     + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt
     3     + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t 
     4     + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt
     5     + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t 
     6     + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt
     7     + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t 
     8     + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt

!..end of statement function definitions

         ht => eos_ht
         
         h = planck_h
         third  = 1.0d0/3.0d0
         sioncon = (2.0d0 * pi * amu * kerg)/(h*h)
         sifac  = 8.6322745944370191d-45
         kergavo = kerg * avo
         asoli3  = asol/3.0d0
         clight2 = clight*clight
         eostol = 1.0d-13
         fpmin  = 1.0d-14
         !..note: sifac = h**3/(2.0d0*pi*amu)**1.5d0
         forth   = 4.0d0/3.0d0
         fiveth  = 5.0d0/3.0d0
         teninth = 10.0d0/9.0d0
         esqu    = qe*qe
         forthpi = forth * pi

         ierr = 0
         
         abar = abar_in
         zbar = zbar_in
         temp = temp_in
         logtemp = logtemp_in
         den = den_in
         logden = logden_in

!..for very low T, convert all H to H2.  adjust abar and zbar accordingly.
         
         ! NOTE: table lookup uses din rather than den
         ytot1 = 1.0d0/abar
         ye    = ytot1 * zbar
         din     = ye*den
         
         skip_elec_pos = must_skip_elec_pos
         if (.not. skip_elec_pos) then ! see if need to set it true
         
            if (temp < ht% templo) then
               if (din > -5d0) then ! clip T so can keep elec_pos
                  temp = ht% templo
                  logtemp = log10(temp)
               else
                  skip_elec_pos = .true.
               end if
            end if
         
            if (din < ht% denlo) then
               skip_elec_pos = .true.
            end if
         
         end if

         if (temp > ht% temphi) then            
            temp = ht% temphi
            logtemp = ht% logthi
         end if
         
         if (din > ht% denhi) then
            din = ht% denhi
         end if
         
         if (skip_elec_pos) then
            abar = 1d0 / (1/abar - Xfrac/2)
            zbar = 1d-10 ! don't set it to 0
            ytot1 = 1.0d0/abar
            ye    = ytot1 * zbar
         end if

!..very neutron rich compositions may need to be bounded, 
!..avoid that extrema for now in order to increase efficiency
c       ye    = max(1.0d-16, ye)


!..initialize local variables
       include 'helm_initialize_local_variables.dek'

!..radiation section:
!       include 'helm_radiation.dek'

!..ion section:
!       include 'helm_ideal_ions.dek'

!..electron-positron section:
       if (.not. skip_elec_pos) then
         include 'helm_electron_positron.dek'
       else ! drop the electron-positron section at very low T
         pele    = 0.0d0
         dpeledd = 0.0d0
         dpeledt = 0.0d0
         dpeleda = 0.0d0
         dpeledz = 0.0d0
         eele    = 0.0d0
         deeledd = 0.0d0
         deeledt = 0.0d0
         deeleda = 0.0d0
         deeledz = 0.0d0
         sele    = 0.0d0
         dseledd = 0.0d0
         dseledt = 0.0d0
         dseleda = 0.0d0
         dseledz = 0.0d0
         ppos    = 0.0d0
         dpposdd = 0.0d0
         dpposdt = 0.0d0
         dpposda = 0.0d0
         dpposdz = 0.0d0
         epos    = 0.0d0
         deposdd = 0.0d0
         deposdt = 0.0d0
         deposda = 0.0d0
         deposdz = 0.0d0
         spos    = 0.0d0
         dsposdd = 0.0d0
         dsposdt = 0.0d0
         dsposda = 0.0d0
         dsposdz = 0.0d0
         etaele = -20d0
       end if

!..coulomb section:
       if ((ht% with_coulomb_corrections) .and. (.not. skip_elec_pos)) then
!         include 'helm_coulomb2.dek'
       else
         pcoul    = 0.0d0
         dpcouldd = 0.0d0
         dpcouldt = 0.0d0
         dpcoulda = 0.0d0
         dpcouldz = 0.0d0
         ecoul    = 0.0d0
         decouldd = 0.0d0
         decouldt = 0.0d0
         decoulda = 0.0d0
         decouldz = 0.0d0
         scoul    = 0.0d0
         dscouldd = 0.0d0
         dscouldt = 0.0d0
         dscoulda = 0.0d0
         dscouldz = 0.0d0
         plasg = 0
       end if

!..sum the gas and total (gas + radiation) components
       include 'helm_sum_totals.dek'

!..compute the derivative quantities (cv, gamma1 ...etc)
!      include 'helm_gammas.dek'

!..maxwell relations; each is zero if the consistency is perfect
!..if you don't need this, save three divides and comment this out
       ! x   = den * den
       ! dse = temp*dentrdt/denerdt - 1.0d0
       ! dpe = (denerdd*x + temp*dpresdt)/pres - 1.0d0
       ! dsp = -dentrdd*x/dpresdt - 1.0d0
       
       if (.false.) then
       !write(*,'(a30,1pe26.16)') 'temp', temp
       !write(*,'(a30,1pe26.16)') 'den', den
       !write(*,'(a30,1pe26.16)') 'dentrdt', dentrdt
       !write(*,'(a30,1pe26.16)') 'denerdt', denerdt
       !write(*,'(a30,1pe26.16)') 'denerdd', denerdd
       !write(*,'(a30,1pe26.16)') 'dpresdt', dpresdt
       !write(*,'(a30,1pe26.16)') 'dentrdd', dentrdd
       !write(*,'(a30,1pe26.16)') 'dsraddd', dsraddd
       !write(*,'(a30,1pe26.16)') 'dsgasdd', dsgasdd
       write(*,'(a30,1pe26.16)') 'dsiondd', dsiondd
       write(*,'(a30,1pe26.16)') 'dsepdd', dsepdd
       write(*,'(a30,1pe26.16)') 'dscouldd', dscouldd
       write(*,'(a30,1pe26.16)') 'plasg', plasg
       
       write(*,*) 'ht% with_coulomb_corrections', ht% with_coulomb_corrections
       write(*,*) 'skip_elec_pos', skip_elec_pos
 
       stop
       end if

!..store results
      include 'helm_store_results.dek'
      helm_res(h_crp) = sion

!..debugging printout      
      if (.false.) then
         include 'helm_print_results.dek'
      end if

      return
      end subroutine helmeos2aux


      end module
      


