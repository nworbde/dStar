program test_constants
	use constants_def
    use constants_lib
    integer :: ierr
            
    call constants_init('',ierr)
    if (ierr /= 0) stop
    
	call do_one('pi',pi)
	call do_one('twopi',twopi)
	call do_one('threepi',threepi)
	call do_one('fourpi',fourpi)
	call do_one('onethird',onethird)
	call do_one('twothird',twothird)
	call do_one('fivethird',fivethird)
	call do_one('seventhird',seventhird)
	call do_one('threepisquare',threepisquare)
	call do_one('boltzmann',boltzmann)
	call do_one('avogadro',avogadro)
	call do_one('clight',clight)
	call do_one('clight2',clight2)
	call do_one('amu',amu)
	call do_one('hbar',hbar)
	call do_one('finestructure',finestructure)
	call do_one('electroncharge',electroncharge)
	call do_one('GMsun',GMsun)
	call do_one('Gnewton',Gnewton)
	call do_one('Msun',Msun)
	call do_one('Mneutron',Mneutron)
	call do_one('Mproton',Mproton)
	call do_one('Melectron',Melectron)
	call do_one('Mmuon',Mmuon)
	call do_one('mev_to_ergs',mev_to_ergs)
	call do_one('fm_to_cm',fm_to_cm)
	call do_one('ergs_to_mev',ergs_to_mev)
	call do_one('cm_to_fm',cm_to_fm)
	call do_one('a_Bohr',a_Bohr)
    call do_one('a_Bohr/a_Bohr(CODATA)-1',a_Bohr/5.29177210903e-9_dp-1.0_dp)
	call do_one('Rydberg',Rydberg)
    call do_one('Rydberg/Rydberg(CODATA)-1',Rydberg/2.1798723611035e-11_dp-1.0_dp)
    call do_one('Thomson',Thomson)
    call do_one('Thomson/Thomson(CODATA)-1',Thomson/6.6524587321e-25_dp-1.0_dp)
	call do_one('hbarc_n',hbarc_n)
	call do_one('mass_n',mass_n)
	call do_one('length_n',length_n)
	call do_one('density_n',density_n)
	call do_one('pressure_n',pressure_n)
	call do_one('temperature_n',temperature_n)
	call do_one('Mn_n',Mn_n)
	call do_one('Mp_n',Mp_n)
	call do_one('Me_n',Me_n)
	call do_one('mass_g',mass_g)
	call do_one('potential_g',potential_g)
	call do_one('length_g',length_g)
	call do_one('density_g',density_g)
	call do_one('pressure_g',pressure_g)
	call do_one('time_g',time_g)
	
	contains
	subroutine do_one(name,val)
		character(len=*), parameter :: fmt = '(a28," = ",es20.13)'
		character(len=*), intent(in) ::  name
		real(dp), intent(in) :: val
		write (*,fmt) name, val
	end subroutine do_one

end program test_constants
