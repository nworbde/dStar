module HZ90_comp
    use math_lib
    use const_def, only: dp
    use exceptions_lib
    use nucchem_def
    use nucchem_lib
    implicit none
    integer, parameter :: HZ90_number_isos = 19
    character(len=iso_name_length), parameter, dimension(HZ90_number_isos) :: HZ90_network = &
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
    integer, parameter :: HZ90_number_transitions = 17
    integer, parameter :: HZ90_number_layers = HZ90_number_transitions+1
    character(len=iso_name_length), parameter, dimension(HZ90_number_layers) :: HZ90_ion_composition = [ &
    &   character(len=iso_name_length) :: &
    &   'fe56','cr56','ti56','ca56','ar56','s52','si46','mg40','ca68','ar62','s56','si50', &
    &   'mg44','ar66','s60','si54','mg48','ti88' ]
    real(dp),parameter, dimension(HZ90_number_transitions) :: HZ90_transition_pressures = [ &
    &   7.235d26,9.569d27,1.152d29,4.747d29,1.361d30,1.980d30,2.253d30, &
    &   2.637d30,2.771d30,3.216d30,3.825d30,4.699d30,6.043d30,7.233d30,9.238d30, &
    &   1.228d31,1.602d31 ]
    real(dp), parameter, dimension(HZ90_number_layers) :: Xn = [ &
    &   0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,  &
    &   4.0_dp/56.0_dp, 10.0_dp/56.0_dp, 16.0_dp/56.0_dp, 22.0_dp/56.0_dp, &
    &   50.0_dp/112.0_dp, 56.0_dp/112.0_dp, 62.0_dp/112.0_dp, 68.0_dp/112.0_dp, 158.0_dp/224.0_dp, &
    &   164.0_dp/224.0_dp, 170.0_dp/224.0_dp, 176.0_dp/224.0_dp, 360.0_dp/448.0_dp ]

contains
    
    subroutine do_generate_HZ90_table(lgP,Y)
        real(dp), dimension(:), intent(in) :: lgP
        real(dp), dimension(HZ90_number_isos,size(lgP)), intent(out) :: Y
    
        integer, dimension(max_nnuclib) :: network_indcs
        integer :: nz, nion, i, indx, n_indx, ierr
        integer, dimension(HZ90_number_isos) :: indcs
        real(dp), dimension(HZ90_number_transitions) :: lgPt
        
        ! local variables
        lgPt = log10(HZ90_transition_pressures)
        nion = HZ90_number_isos
    
        ! indcs are locations of nuclei in nuclib database
        indcs = [(get_nuclide_index(HZ90_network(i)),i=1,HZ90_number_isos)]
        ! network_indcs allow for reverse lookup of specific nuclei
        network_indcs = 0
        network_indcs(indcs) = [(i,i=1,HZ90_number_isos)]
        n_indx = network_indcs(get_nuclide_index('n'))
        
        ! loop over layers to set composition
        Y = 0.0_dp
        ! first layer
        indx = network_indcs(get_nuclide_index(HZ90_ion_composition(1)))
        where(lgP <= lgPt(1)) 
            Y(n_indx,:) = Xn(1)
            Y(indx,:) = (1.0_dp-Xn(1))/nuclib% A(indcs(indx))
        end where
        ! middle layers
        do i = 2, HZ90_number_transitions
            indx = network_indcs(get_nuclide_index(HZ90_ion_composition(i)))
            where(lgP > lgPt(i-1) .and. lgP <= lgPt(i))
                Y(n_indx,:) = Xn(i)
                Y(indx,:) = (1.0_dp-Xn(i))/nuclib% A(indcs(indx))
            end where
        end do
        ! last layer
        indx = network_indcs(get_nuclide_index(HZ90_ion_composition(HZ90_number_layers)))
        where (lgP > lgPt(HZ90_number_transitions))
            Y(n_indx,:) = Xn(HZ90_number_layers)
            Y(indx,:) = (1.0_dp-Xn(HZ90_number_layers))/nuclib% A(indcs(indx))
        end where
        
    end subroutine do_generate_HZ90_table

end module HZ90_comp
