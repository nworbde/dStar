module nucchem_def
    use constants_def, only: dp
    use utils_def, only: integer_dict

    implicit none

    integer, parameter :: nuclide_not_found = -1
    ! max. Z in the nucchem database
    integer, parameter :: max_element_z = 85
    
    ! data storage parameter for nuclib
    ! each isotope has a name and "provenance" -- a reference to where the 
    ! mass, spin, partition function are described.
    integer, parameter :: iso_name_length = 5
    integer, parameter :: provenance_length = 6
    integer, parameter :: reaction_reference_length = 4
    ! maximum number of nuclides in nucchem database
    integer, parameter :: max_nnuclib=10000
    ! no. entries in partition fcn table
    integer, parameter :: npfcn = 24
    
    ! for storage of composition information
    type composition_info_type
        real(dp) :: A
        real(dp) :: Z
        real(dp) :: Z53
        real(dp) :: Z2
        real(dp) :: Z73
        real(dp) :: Z52
        real(dp) :: ZZ1_32
        real(dp) :: Z2XoA2
        real(dp) :: Ye
        real(dp) :: Yn
        real(dp) :: Q
    end type composition_info_type
    
    ! table of elements; note that hydrogen is referred to as 'p' in the table 
    ! and neutrons are referred to as 'n'
    character(len=iso_name_length), dimension(0:max_element_z) ::  &
    & element_name = [character(len=iso_name_length) :: & 
    & 'neut','h','he','li','be','b','c','n','o','f','ne', 'na','mg','al','si', &
    & 'p', 's', 'cl', 'ar', 'k', 'ca', 'sc', 'ti', 'v','cr','mn', 'fe', 'co', &
    & 'ni','cu','zn','ga','ge','as','se','br','kr','rb','sr','y','zr', &
    & 'nb','mo','tc','ru','rh','pd','ag','cd','in','sn','sb','te','i','xe', &
    & 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy','ho','er', &
    & 'tm','yb', 'lu','hf','ta','w','re','os','ir','pt','au','hg', &
    & 'tl','pb','bi','po','at' ]
    
    ! synonyms of hydrogen isotopes, index by mass number
    character(len=iso_name_length), dimension(3) ::  &
    & h_isotopes = [character(len=iso_name_length) :: 'p','d','t']

    ! aluminum-26 isomers
    character(len=iso_name_length), dimension(2:3) ::  &
    & al_isomers = [character(len=iso_name_length) :: 'al-6','al*6']
    
    ! storage for nuclib database
    type nuclib_data
        integer :: Nnuclides
        ! for each nuclide, store its name, provenance, mass (real), charge and 
        ! neutron numbers, ground-state spin, mass excess (MeV) and partition 
        ! function table
        character(len=iso_name_length), dimension(:), allocatable :: name
        character(len=provenance_length),dimension(:), allocatable :: provenance
        real(dp), dimension(:), allocatable :: A
        integer, dimension(:), allocatable  :: Z
        integer, dimension(:), allocatable :: N
        real(dp), dimension(:), allocatable :: spin
        real(dp), dimension(:), allocatable :: mass_excess
        real(dp), dimension(:,:), allocatable :: pfcn
    end type nuclib_data
    
    ! temperatures, partition function evaluation
    real(dp), dimension(npfcn), save :: pfcn_T9

    type(integer_dict), pointer, save :: nuclide_dict=>null()
    type(nuclib_data), save :: nuclib
    logical, save :: nucchem_is_initialized
    
end module nucchem_def
