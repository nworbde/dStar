module composition_models

    use const_def, only: dp
    
    real(dp), parameter :: transition_width = 0.02  
    
    integer, parameter :: HZ90_number = 19
    character(len=iso_name_length), parameter, dimension(HZ90_number) :: HZ90_network = [ character(len=iso_name_length) :: &
    &   'n', &
    &   'mg40', &
    &   'mg44', &
    &   'mg48', &
    &   'si46', &
    &   'si50' &
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
    &   'ti88'  &
    &   'cr56', &
    &   'fe56' ]
    
contains
    
    subroutine HZ90(lgrho, X)
        use nucchem_def
        real(dp), intent(in), dimension(:) :: lgrho
        real(dp), intent(out), dimension(:,:), allocatable :: X
        integer, dimension(max_nnuclib) :: network_indcs
         
        integer :: Ntab, i, indx
        integer, dimension(HZ90_number :: indcs
        integer, parameter :: number_layers = 17
        character(len=iso_name_length), parameter, dimension(number_layers+1) :: ion_composition = [ &
        &   character(len=iso_name_length) :: &
        &   'fe56','cr56','ti56','ca56','ar56','s52','si46','mg40','ca68','ar62','s56','si50', &
        &   'mg44','ar66','s60','si54','mg48','ti88' ]
        
        real(dp),parameter, dimension(number_layers) :: transition_densities = [ &
        &   1.494d9,1.114d10,7.848d10,2.496d11,6.110d11,9.075d11,1.131d12,1.455d12,1.766d12,2.134d12, &
        &   2.634d12,3.338d124.379d12,5.839d12,7.041d12,8.980d12,1.127d13 ]

        real(dp), dimension(number_layers) :: lg_rhot
        real(dp), parameter, dimension(number_layers) :: Xn = [ &
        &   0.0, 0.0, 0.0, 0.0, 0.0, 0.07, 0.18, 0.29, 0.39, 0.45, 0.50, 0.55, 0.61, 0.70, 0.73, 0.76, 0.80]
        
        Ntab = length(lgrho)
        allocate(lgP(Ntab), X(HZ90_number,Ntab))
        
        lg_rhot = log10(transition_densities)
        
        ! set the network pointers
        indcs = [(get_nuclide_index(HZ90_network(i)),i=1,HZ90_number)]

        ! and the reverse lookup
        network_indcs = 0
        network_indcs(indcs) = [(i,i=1,HZ90_number)]
        
        ! for each layer set composiiton...we'll smooth the transitions in the next step
        X = 0.0
        ! first layer
        indx = network_indcs(get_nuclide_index(ion_composition(1)))
        where(lgrho <= lg_rhot(1)) 
            X(indx,:) = 1.0
        end where
        
        do i = 2, number_layers
            indx = network_indcs(get_nuclide_index(ion_composition(i)))
            where(lgrho > lg_rhot(i-1) .and. lgrho <= lg_rhot(i))
                X(indx,:) = 1.0
            end where
        end do
        
        indx = network_indcs(get_nuclide_index(ion_composition(number_layers+1)))
        where (lgrho > lg_rhot(number_layers))
            X(indx,:) = 1.0
        end where
        
        ! now smooth the transitions
        do i = 1, number_layers
            lgrho1 = lg_rhot - transition_width
            lgrho2 = lg_rhot + transition_width
            width = lgrho2-lgrho1
            indx1 = network_indcs(get_nuclide_index(ion_composition(i)))
            indx2 = network_indcs(get_nuclide_index(ion_composition(i+1)))
            where(lgrho >= lgrho1 .and. lgrho <= lgrho2)
                X(indx1,:) = cos(0.5*pi*(lgrho-lgrho1/width))
                X(indx2,:) = 1.0 - cos(0.5*pi*(lgrho-lgrho1/width))
            end where
        end do
        
    end subroutine HZ90
    

end module composition_models
