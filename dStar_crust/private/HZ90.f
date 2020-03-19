module hz90

    use const_def, only: dp
    use nucchem_def, only: iso_name_length
    
    integer, parameter :: HZ90_number = 19
    character(len=iso_name_length), parameter, dimension(HZ90_number) :: HZ90_network = &
    &   [ character(len=iso_name_length) :: &
    &   'n', &
    &   'mg40', &
    &   'mg44', &
    &   'mg48', &
    &   'si46', &
    &   'si50', &
    &   'si54', &
    &   's52', &
    &   's56', &
    &   's60', &
    &   'ar56', &
    &   'ar62', &
    &   'ar66', &
    &   'ca56', &
    &   'ca68', &
    &   'ti56', &
    &   'ti88',  &
    &   'cr56', &
    &   'fe56' ]
    integer, parameter :: number_layers = 17
    character(len=iso_name_length), parameter, dimension(number_layers+1) :: ion_composition = [ &
    &   character(len=iso_name_length) :: &
    &   'fe56','cr56','ti56','ca56','ar56','s52','si46','mg40','ca68','ar62','s56','si50', &
    &   'mg44','ar66','s60','si54','mg48','ti88' ]
    real(dp),parameter, dimension(number_layers) :: transition_pressures = [ &
    &   7.235d26,9.569d27,1.152d29,4.747d29,1.361d30,1.980d30,2.253d30, &
    &   2.637d30,2.771d30,3.216d30,3.825d30,4.699d30,6.043d30,7.233d30,9.238d30, &
    &   1.228d31,1.602d31 ]
    real(dp), parameter, dimension(number_layers+1) :: Xn = [ &
    &   0.0, 0.0, 0.0, 0.0, 0.0,  &
    &   4.0/56.0, 10.0/56.0, 16.0/56.0, 22.0/56.0, &
    &   50.0/112.0, 56.0/112.0, 62.0/112.0, 68.0/112.0, 158.0/224.0, &
    &   164.0/224.0, 170.0/224.0, 176.0/224.0, 360.0/448.0 ]
    
    integer, parameter :: HZ08_number = 34
    character(len=iso_name_length), parameter, dimension(HZ08_number) :: HZ08_network =  &
    &   [ character(len=iso_name_length) :: &
    &   'n', &
    &   'ne36', &
    &   'mg46', &
    &   'mg42', &
    &   'si50', &
    &   'si54', &
    &   'si62', &
    &   's56', &
    &   's60', &
    &   's68', &
    &   'ar62', &
    &   'ar66', &
    &   'ar74', &
    &   'ca68', &
    &   'ca72', &
    &   'ca80', &
    &   'ti74', &
    &   'ti86', &
    &   'ti116', &
    &   'cr80', &
    &   'cr92', &
    &   'cr118', &
    &   'fe86', &
    &   'fe120', &
    &   'ni92', &
    &   'ni124', &
    &   'ge106', &
    &   'se106', &
    &   'kr106', &
    &   'sr106', &
    &   'zr106', &
    &   'mo106', &
    &   'ru106', &
    &   'pd106' ]
    integer, parameter :: number_HZ08_layers = 29    
    character(len=iso_name_length), parameter, dimension(number_HZ08_layers+1) :: HZ08_ion_composition = [ &
    &   character(len=iso_name_length) ::  &
    &   'pd106', 'ru106', 'mo106', 'zr106', 'sr106', 'kr106', 'se106', 'ge106', 'ni92',  &
    &   'fe86', 'cr80', 'ti74','ca68', 'ar62', 's56', 'si50', 'mg42', 'ca72', 'ar66', 's60', &
    &   'si54', 'cr92', 'ti86', 'ca80', 'ar74', 's68', 'ni124', 'fe120', 'cr118', 'ti116' ]
!     real(dp), parameter, dimension(number_HZ08_layers) :: HZ08_transition_pressures = [ &
!     & ]
!     real(dp), parameter, dimension(number_HZ08_layers+1) :: HZ08_Xn = [ &
!     & ]
contains
    
    subroutine do_make_crust(lgP, Yion, Xneut, charged_ids, ncharged, ion_info)
        use const_def
        use nucchem_def
        use nucchem_lib
        real(dp), intent(in), dimension(:) :: lgP
        real(dp), intent(out), dimension(:,:) :: Yion   ! (HZ90_number, size(lgP))
        real(dp), intent(out), dimension(:) :: Xneut    ! (size(lgP))
        integer, intent(out), dimension(HZ90_number) :: charged_ids
        integer, intent(out) :: ncharged
        type(composition_info_type), dimension(:), intent(out) :: ion_info   ! (size(lgP))
        integer, dimension(max_nnuclib) :: network_indcs
        
        real(dp), allocatable, dimension(:,:) :: X
        integer :: Ntab, i, indx, n_indx, indx1, indx2
        integer, dimension(HZ90_number) :: indcs
        
        real(dp), dimension(number_layers) :: lg_Pt        
        real(dp) :: lgP1, lgP2, width, Xsum
        
        Ntab = size(lgP)
        allocate(X(HZ90_number,Ntab))
        
        lg_Pt = log10(transition_pressures)
        
        ! set the network pointers
        indcs = [(get_nuclide_index(HZ90_network(i)),i=1,HZ90_number)]

        ! and the reverse lookup
        network_indcs = 0
        network_indcs(indcs) = [(i,i=1,HZ90_number)]
        n_indx = network_indcs(get_nuclide_index('n'))
        
        ! for each layer set composiiton
        X = 0.0
        ! first layer
        indx = network_indcs(get_nuclide_index(ion_composition(1)))
        where(lgP <= lg_Pt(1)) 
            X(n_indx,:) = Xn(1)
            X(indx,:) = 1.0-Xn(1)
        end where
        ! middle layers
        do i = 2, number_layers
            indx = network_indcs(get_nuclide_index(ion_composition(i)))
            where(lgP > lg_Pt(i-1) .and. lgP <= lg_Pt(i))
                X(n_indx,:) = Xn(i)
                X(indx,:) = 1.0-Xn(i)
            end where
        end do
        ! last layer
        indx = network_indcs(get_nuclide_index(ion_composition(number_layers+1)))
        where (lgP > lg_Pt(number_layers))
            X(n_indx,:) = Xn(number_layers+1)
            X(indx,:) = 1.0-Xn(number_layers+1)
        end where
        
        ! composition moments
        do i = 1, Ntab
            call compute_composition_moments(HZ90_number,indcs,X(:,i), &
            &   ion_info(i),Xsum,ncharged,charged_ids,Yion(:,i),exclude_neutrons=.TRUE., &
            &   abunds_are_mass_fractions=.TRUE.)
        end do
        Xneut = X(n_indx,:)
        
        deallocate(X)
    end subroutine do_make_crust
    
    subroutine find_densities(eos_handle,lgP,lgRho,lgEps,Yion,ncharged,charged_ids,ionic,Tref)
        use exceptions_lib
        use constants_def
        use nucchem_def
        use num_lib

        real(dp) :: Pfac
        integer, intent(in) :: eos_handle
        real(dp), dimension(:), intent(in) :: lgP
        real(dp), dimension(:), intent(out) :: lgRho
        real(dp), dimension(:), intent(out) :: lgEps
        real(dp), dimension(:,:), intent(in) :: Yion
        integer, intent(in) :: ncharged
        integer, dimension(:), intent(in) :: charged_ids
        real(dp), intent(in) :: Tref
        type(composition_info_type), dimension(:), intent(in) :: ionic
        real(dp), dimension(:), pointer :: rpar=>null()
        integer, dimension(:), pointer :: ipar=>null()
        integer :: lipar, lrpar
        integer :: i,Ntab
        real(dp) :: x1, x3, y1, y3, epsx, epsy, lgRho_guess
        integer :: imax, ierr
        type(warning) :: rootfind_error=warning(scope='find_densities')
        type(alert) :: rootfind_status=alert(scope='find_densities')
        character(len=128) :: error_message
        
        Pfac = 0.25*(threepisquare)**onethird *hbar*clight*avo**(4.0*onethird)
        Ntab = size(lgP)
        imax = 20
        epsx = 1.0d-12
        epsy = 1.0d-12
        
        ! decide the size of the parameter arrays
        lipar = 2 + ncharged
        allocate(ipar(lipar))
        ipar(1) = eos_handle
        ipar(2) = ncharged
        ipar(3:ncharged+2) = charged_ids(1:ncharged)
        
        lrpar = ncharged + 11 + 4
        allocate(rpar(lrpar))
        
        ! last value of rpar is a guess for the density; if 0, will be calculated for relativistic electron gas
        rpar(lrpar) = 0.0
        do i = 1, Ntab
            ! stuff the composition information into the parameter array
            rpar(1:ncharged) = Yion(1:ncharged,i)
            rpar(ncharged+1) = ionic(i)% A
            rpar(ncharged+2) = ionic(i)% Z
            rpar(ncharged+3) = ionic(i)% Z53
            rpar(ncharged+4) = ionic(i)% Z2
            rpar(ncharged+5) = ionic(i)% Z73
            rpar(ncharged+6) = ionic(i)% Z52
            rpar(ncharged+7) = ionic(i)% ZZ1_32
            rpar(ncharged+8) = ionic(i)% Z2XoA2
            rpar(ncharged+9) = ionic(i)% Ye
            rpar(ncharged+10) = ionic(i)% Yn
            rpar(ncharged+11) = ionic(i)% Q
            rpar(ncharged+12) = lgP(i)
            rpar(ncharged+13) = Tref
            
            if (i > 1 .and. rpar(lrpar) /= 0.0) then
                lgRho_guess = lgRho(i-1) + (lgP(i)-lgP(i-1))/rpar(lrpar)
            else
                lgRho_guess = log10((10.0**lgP(i)/Pfac)**0.75 /ionic(i)% Ye)
            end if
            
            call look_for_brackets(lgRho_guess,0.05*lgRho_guess,x1,x3,match_density, &
            &   y1,y3,imax,lrpar,rpar,lipar,ipar,ierr)
            if (rootfind_error% raised(ierr)) then
                write(error_message,'(a,5f0.3)') &
                &   'unable to bracket root: ',lgP(i), x1, x3, y1, y3
                call rootfind_status% report(error_message)
                cycle
            end if
            
            lgRho(i) = safe_root_with_initial_guess(match_density,lgRho_guess,x1,x3,y1,y3, &
            &   imax,epsx,epsy,lrpar,rpar,lipar,ipar,ierr)
            if (rootfind_error% raised(ierr)) then
                write(error_message,'(a,5f0.3)') &
                &   'unable to converge: ',lgP(i), x1, x3, y1, y3
                call rootfind_status% report(error_message)
                cycle
            end if
            lgEps(i) = rpar(ncharged+14)
        end do
    end subroutine find_densities
    
    real(dp) function match_density(lgRho, dfdlgRho, lrpar, rpar, lipar, ipar, ierr)
       ! returns with ierr = 0 if was able to evaluate f and df/dx at x
       ! if df/dx not available, it is okay to set it to 0
       use ieee_arithmetic
       use exceptions_lib
       use constants_def
       use superfluid_def, only: max_number_sf_types
       use superfluid_lib, only: sf_get_results
       use nucchem_def
       use dStar_eos_def
       use dStar_eos_lib
       
       integer, intent(in) :: lrpar, lipar
       real(dp), intent(in) :: lgRho
       real(dp), intent(out) :: dfdlgRho
       integer, intent(inout), pointer :: ipar(:) ! (lipar)
       real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
       integer, intent(out) :: ierr
       integer :: eos_handle, ncharged
       type(composition_info_type) :: ionic
       integer, dimension(:), allocatable :: charged_ids
       real(dp), dimension(:), allocatable :: Yion
       real(dp), dimension(num_dStar_eos_results) :: res
       integer :: phase
       real(dp) :: chi, lgPwant, lgP, kFn, kFp, Tcs(max_number_sf_types)
       real(dp) :: rho, T, Eint
       type(assertion) :: lgP_is_num = assertion(scope='match_density', &
       &    message='pressure is a number')
       type(crust_eos_component), dimension(num_crust_eos_components) :: eos_components
    
       ierr = 0       
       eos_handle = ipar(1)
       ncharged = ipar(2)
       allocate(charged_ids(ncharged),Yion(ncharged))
       charged_ids(:) = ipar(3:ncharged+2)
       
       Yion = rpar(1:ncharged)
       ionic% A = rpar(ncharged+1)
       ionic% Z = rpar(ncharged+2)
       ionic% Z53 = rpar(ncharged+3)
       ionic% Z2 = rpar(ncharged+4)
       ionic% Z73 = rpar(ncharged+5)
       ionic% Z52 = rpar(ncharged+6)
       ionic% ZZ1_32 = rpar(ncharged+7)
       ionic% Z2XoA2 = rpar(ncharged+8)
       ionic% Ye = rpar(ncharged+9)
       ionic% Yn = rpar(ncharged+10)
       ionic% Q = rpar(ncharged+11)
       
       rho = 10.0**lgRho
       T = rpar(ncharged+13)
       ! when searching for a root, the density can be > nuclear, so that chi > 1.0,
       ! which results in a bad neutron number density.
       ! we cap chi at the value appropriate a density of 1.0e14 g/cc.
       chi = nuclear_volume_fraction(min(rho,1.0e14_dp),ionic,default_nuclear_radius)

       kFp = 0.0_dp
       kFn = neutron_wavenumber(rho,ionic,chi)
       call sf_get_results(kFp,kFn,Tcs)
       call eval_crust_eos(eos_handle,rho,T,ionic,ncharged, &
       &    charged_ids,Yion,Tcs,res,phase,chi,eos_components)
       Eint = res(i_lnE)
       
       lgPwant = rpar(ncharged+12)
       lgP = res(i_lnP)/ln10
       call lgP_is_num% assert(.not. ieee_is_nan(lgP))
       
       ! mass-energy density
       rpar(ncharged+14) = log10(rho* &
       &    (1.0+dot_product(Yion(:),nuclib%mass_excess(charged_ids))/amu_n + Eint/clight2))
       rpar(lrpar) = res(i_chiRho)
       match_density = lgP - lgPwant
       deallocate(charged_ids,Yion)
    end function match_density
    

end module hz90
