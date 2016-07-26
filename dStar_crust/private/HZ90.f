module hz90

    use const_def, only: dp
    use nucchem_def, only: iso_name_length
    
    real(dp), parameter :: transition_width = 0.02  
    
    integer, parameter :: HZ90_number = 19
    integer, parameter :: HZ08_number = 34
    character(len=iso_name_length), parameter, dimension(HZ90_number) :: HZ90_network = [ character(len=iso_name_length) :: &
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
    
    character(len=iso_name_length), parameter, dimension(HZ08_number) :: HZ08_network = [ character(len=iso_name_length) :: &
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
    
    integer, parameter :: number_layers = 17
    
    character(len=iso_name_length), parameter, dimension(number_layers+1) :: ion_composition = [ &
    &   character(len=iso_name_length) :: &
    &   'fe56','cr56','ti56','ca56','ar56','s52','si46','mg40','ca68','ar62','s56','si50', &
    &   'mg44','ar66','s60','si54','mg48','ti88' ]
    
    real(dp),parameter, dimension(number_layers) :: transition_pressures = [ &
    &   7.235d26,9.569d27,1.152d29,4.747d29,1.361d30,1.980d30,2.253d30,2.637d30,2.771d30,3.216d30, &
    &   3.825d30,4.699d30,6.043d30,7.233d30,9.238d30,1.228d31,1.602d31 ]

    real(dp), parameter, dimension(number_layers+1) :: Xn = [ &
    &   0.0, 0.0, 0.0, 0.0, 0.0,  &
    &   0.07, 0.18, 0.29, 0.39, 0.45, 0.50, 0.55, 0.61, 0.70, 0.73, 0.76, 0.80, 0.80]

    
    integer, parameter :: number_HZ08_layers = 29
    
    character(len=iso_name_length), parameter, dimension(number_HZ08_layers+1) :: HZ08_ion_composition = [ &
    &   character(len=iso_name_length) ::  &
    &   'pd106', 'ru106', 'mo106', 'zr106', 'sr106', 'kr106', 'se106', 'ge106', 'ni92', 'fe86', 'cr80', 'ti74', &
    &   'ca68', 'ar62', 's56', 'si50', 'mg42', 'ca72', 'ar66', 's60', 'si54', 'cr92', 'ti86', 'ca80', 'ar74', &
    &   's68', 'ni124', 'fe120', 'cr118', 'ti116' ]
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
        
        ! for each layer set composiiton...we'll smooth the transitions in the next step
        X = 0.0
        ! first layer
        indx = network_indcs(get_nuclide_index(ion_composition(1)))
        n_indx = network_indcs(get_nuclide_index('n'))
        where(lgP <= lg_Pt(1)) 
            X(n_indx,:) = Xn(1)
            X(indx,:) = 1.0-Xn(1)
        end where
        
        do i = 2, number_layers
            indx = network_indcs(get_nuclide_index(ion_composition(i)))
            where(lgP > lg_Pt(i-1) .and. lgP <= lg_Pt(i))
                X(n_indx,:) = Xn(i)
                X(indx,:) = 1.0-Xn(i)
            end where
        end do
        
        indx = network_indcs(get_nuclide_index(ion_composition(number_layers+1)))
        where (lgP > lg_Pt(number_layers))
            X(n_indx,:) = Xn(number_layers+1)
            X(indx,:) = 1.0-Xn(number_layers+1)
        end where
        
        ! now smooth the transitions
        do i = 1, number_layers
            lgP1 = lg_Pt(i) - transition_width
            lgP2 = lg_Pt(i) + transition_width
            width = 2.0*transition_width
            indx1 = network_indcs(get_nuclide_index(ion_composition(i)))
            indx2 = network_indcs(get_nuclide_index(ion_composition(i+1)))
            where(lgP >= lgP1 .and. lgP <= lgP2)
                X(n_indx,:) = (Xn(i)-Xn(i+1))*cos(0.5*pi*(lgP-lgP1)/width) + Xn(i+1)
                X(indx1,:) = (1.0-X(n_indx,:))*cos(0.5*pi*(lgP-lgP1)/width)
                X(indx2,:) = (1.0-X(n_indx,:))*(1.0 - cos(0.5*pi*(lgP-lgP1)/width))
            end where
        end do
        
        ! loop over and compute composition moments
        do i = 1, Ntab
            call compute_composition_moments(HZ90_number,indcs,X(:,i), &
            &   ion_info(i),Xsum,ncharged,charged_ids,Yion(:,i),exclude_neutrons=.TRUE., &
            &   abunds_are_mass_fractions=.TRUE.)
        end do
        Xneut = X(n_indx,:)
        
        deallocate(X)
    end subroutine do_make_crust   

end module hz90
