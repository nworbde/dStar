program tabulate_mixture
    use math_lib
    use iso_fortran_env, only: output_unit, error_unit
    use utils_lib, only: StrLowCase
    use exceptions_lib
    use constants_lib
    use nucchem_def
    use nucchem_lib
	use superfluid_def, only: max_number_sf_types, neutron_1S0
	use superfluid_lib
    use dStar_eos_lib
    use conductivity_lib
    use PPP_electron
    
    implicit none
    ! default table limits and resolution in lg[rho/(g cm**-3)]
    integer :: Nlgrho
    real(dp), parameter :: default_lgrho_min=5.0_dp, default_lgrho_max = 9.7_dp, default_delta_lgrho = 0.1_dp
    ! set the table limits and resolution in lg[T/K]
    integer :: NlgT
    real(dp), parameter :: default_lgT_min=6.5_dp, default_lgT_max = 9.0_dp, default_delta_lgT = 0.1_dp
    ! storage for lg(rho), lg(T)
    real(dp), dimension(:), allocatable :: lgrho, lgT
    ! storage for lg(K) (conductivity, both fml and table), eta_e (eletron degeneracy),
    ! Gamma (Coulomb pot. between ions/kT) 
    real(dp), dimension(:,:), allocatable :: lgK_fml, lgK_tab, eta_e, Gamma

    ! user inputs
    ! directory for EOS and composition datafiles
    character(len=128) :: datadir
    ! inlist filename
    character(len=128) :: inlist_filename
    real(dp) :: lgrho_min, lgrho_max, delta_lgrho
    real(dp) :: lgT_min, lgT_max, delta_lgT
    ! composition
    integer, parameter :: max_number_nuclides = 5
    character(len=iso_name_length), dimension(max_number_nuclides) :: nuclides
    real(dp), dimension(max_number_nuclides) :: mass_fractions
    character(len=16) :: method ! 'mean Z', 'rms Z', 'Z2/Z', 'density averaged'
    logical :: add_impurity
    
    ! for interfacing with eos, chem, and conductivity modules
    ! for the nuclide database
    integer :: nnuclides
    integer, dimension(:), allocatable :: chem_ids, charged_ids
    real(dp), dimension(:), allocatable :: X, Y, Yion, Z, A
    integer :: ncharged
    type(composition_info_type) :: ionic
    real(dp) :: chi, Xsum
    ! for the equation of state
    integer :: eos_handle, phase
    real(dp), dimension(num_dStar_eos_results) :: res
    type(crust_eos_component), dimension(num_crust_eos_components) :: eos_components
    ! superfluid critical temperatures
    real(dp), dimension(max_number_sf_types) :: Tcs

    ! for the thermal conductivity
    integer :: cond_handle
    type(conductivity_components) :: kappa

    ! for error checking
    integer :: ierr
    type(assertion) :: check_okay=assertion(scope='main')
    type(assertion) :: open_namelist=assertion(scope='read_controls')
    type(assertion) :: read_namelist=assertion(scope='read_controls')
    type(assertion) :: chemistry_okay=assertion(scope='chemistry')
    type(assertion) :: eos_okay=assertion(scope='eos')
    type(assertion) :: conductivity_okay=assertion(scope='conductivity')
        
    character(len=64) :: progname
    ! miscellaneous
    integer :: i,j,k
    real(dp) :: rho,T,K_e,eta,mu_e
    integer :: iounit

    namelist /controls/ &
        datadir, &
        lgrho_min, &
        lgrho_max, &
        delta_lgrho, &
        lgT_min, &
        lgT_max, &
        delta_lgT, &
        nnuclides, &
        nuclides, &
        mass_fractions, &
        method, &
        add_impurity
            
    ! set defaults
    datadir = '../../data'
    lgrho_min = default_lgrho_min
    lgrho_max = default_lgrho_max
    delta_lgrho = default_delta_lgrho
    lgT_min = default_lgT_min
    lgT_max = default_lgT_max
    delta_lgT = default_delta_lgT
    nnuclides = 2
    nuclides = ''
    nuclides(1) = 'c12'
    nuclides(2) = 'fe56'
    mass_fractions = 0.0_dp
    mass_fractions(1:2) = 0.5_dp
    method = 'mean Z'
    add_impurity = .FALSE.

    ierr = 0
    ! get user inputs
    inlist_filename = 'inlist'
    call open_namelist% set_message('failed to open '//trim(inlist_filename))
    call read_namelist% set_message('failed while reading '//trim(inlist_filename))
    
    if (command_argument_count() == 1) then
        call get_command_argument(1,inlist_filename)
    end if
    
    open (newunit=iounit,file=trim(inlist_filename),  &
    & action='read', delim='quote', status='old', iostat=ierr)
    call open_namelist% assert(ierr==0)
    read(iounit,nml=controls,iostat=ierr)
    close(iounit)
    call read_namelist% assert(ierr == 0)

    ! turn off all status messages except for warnings and failures
    call set_verbosity(1)
    
    ! initialization
    ! physical constants
    call math_init()
    call constants_init('',ierr)
    call check_okay% assert(ierr==0)
    
    ! isotope data
    call nucchem_init(datadir,ierr)
    call chemistry_okay% assert(ierr==0)
    ! set the chemistry
    call chemistry_okay% set_message('locating nuclide information')
    allocate(chem_ids(nnuclides),charged_ids(nnuclides))
    allocate(X(nnuclides),Y(nnuclides),Yion(nnuclides),Z(nnuclides),A(nnuclides))
    do i = 1, nnuclides
        chem_ids(i) = get_nuclide_index(StrLowCase(nuclides(i)))
        call chemistry_okay% assert(chem_ids(i) /= nuclide_not_found)
        X(i) = mass_fractions(i)
    end do
    call compute_composition_moments(nnuclides,chem_ids,X,ionic,Xsum, &
        & ncharged, charged_ids, Yion, exclude_neutrons=.TRUE.,abunds_are_mass_fractions=.TRUE., &
        & renormalize_mass_fractions=.TRUE.)
    Y = Yion
    Z = [ (nuclib% Z(charged_ids(j)), j = 1, ncharged) ]
    A = [ (nuclib% A(charged_ids(j)), j = 1, ncharged) ]
    
    ! equation of state
    call dStar_eos_startup(datadir)
    eos_handle = alloc_dStar_eos_handle(ierr)
    call eos_okay% assert(ierr==0)
    ! suppress the warnings for degenerate ions
    call dStar_eos_set_controls(eos_handle,suppress_warnings=.TRUE.)
    
    ! thermal conductivity
    call conductivity_startup(datadir)
    cond_handle = alloc_conductivity_handle(ierr)
    call conductivity_okay% assert(ierr==0)
    ! set electron-electron scattering formula to Potekhin, Chabrier and Yakovlev (1997)
    call conductivity_set_controls(cond_handle,which_ee_scattering=icond_pcy)

    ! set superfluid temps to a uniform value; we aren't using these, so the precise value is unimportant
    Tcs = 1.0e9_dp
    
    ! define the grid axes
    Nlgrho = nint((lgrho_max-lgrho_min)/delta_lgrho) + 1
    NlgT = nint((lgT_max-lgT_min)/delta_lgT) + 1
    allocate(lgrho(Nlgrho),lgT(NlgT))
    lgrho = linspace(Nlgrho,lgrho_min,lgrho_max)
    lgT = linspace(NlgT,lgT_min,lgT_max)

    ! allocate tables
    allocate(lgK_fml(NlgT,Nlgrho),lgK_tab(NlgT,Nlgrho))
    allocate(eta_e(NlgT,Nlgrho),Gamma(NlgT,Nlgrho))
    
    ! main loop
    density: do i = 1, Nlgrho
        rho = 10.0_dp**lgrho(i)
        temperature: do j = 1, NlgT
            T = 10.0_dp**lgT(j)
            ! evaluate the equation of state
            chi = use_default_nuclear_size
            call eval_crust_eos(eos_handle,rho,T,ionic, &
                & ncharged, charged_ids, Yion, Tcs,  &
                & res, phase, chi, eos_components)
            eta = res(i_Theta) !1.0/TpT
            Gamma(j,i) = res(i_Gamma)
            mu_e = res(i_mu_e)
            eta_e(j,i) = mu_e*mev_to_ergs/boltzmann/T
            call get_thermal_conductivity(cond_handle,rho,T, &
            &   chi,Gamma(j,i),eta,mu_e,ionic,Tcs(neutron_1S0),kappa)
            lgK_fml(j,i) = log10(kappa% electron_total)
            K_e = get_tabulated_conductivity(trim(method),rho,T,ncharged,Y,Z,A,ionic)
            if (kappa%eQ > 0.0_dp .and. add_impurity) K_e = 1.0/(1.0/K_e + 1.0/kappa%eQ)
            lgK_tab(j,i) = log10(K_e)
        end do temperature
    end do density

    ! write out tables for lg(K), mu_e, eta, Gamma
    call write_table('lgK_fml',lgK_fml)
    call write_table('eta_e',eta_e)
    call write_table('Gamma',Gamma)
    call write_table('lgK_tab',lgK_tab)

    ! shutdown
    deallocate(chem_ids,charged_ids,X,Y,Yion,Z,A)
    call conductivity_shutdown
    call dStar_eos_shutdown
    call nucchem_shutdown   

contains
    function get_tabulated_conductivity(method,rho,T,nnuclides,Y,Z,A,ionic) result(K_e)
        character(len=*), intent(in) :: method
        integer, intent(in) :: nnuclides
        real(dp), intent(in) :: rho,T,Y(nnuclides),Z(nnuclides),A(nnuclides)
        type(composition_info_type), intent(in) :: ionic
        integer :: ierr
        real(dp) :: K_e
        real(dp), dimension(:), allocatable :: K_e_species, rho_species
        integer :: i
        type(assertion) :: interp_okay=assertion(scope='get_tabulated_conductivity')

        select case(method)
        case ('mean Z')
            call eval_PPP_electron_table(rho,T,ionic% Z,K_e,ierr)
            call interp_okay% assert(ierr == 0)
        case ('rms Z')
            call eval_PPP_electron_table(rho,T,sqrt(ionic% Z2),K_e,ierr)
            call interp_okay% assert(ierr == 0)
        case ('Z2/Z')
            call eval_PPP_electron_table(rho,T,ionic% Z2/ionic% Z,K_e,ierr)
            call interp_okay% assert(ierr == 0)
        case ('density averaged')
            ! Average over the electron conductivity, assuming that ei scattering is 
            ! proportional to n_i = Y_i*rho/m_u = density of each nuclide. The density 
            ! is adjusted so that n_e is the same for each call.
            allocate(K_e_species(nnuclides),rho_species(nnuclides))
            do i = 1, nnuclides
                rho_species(i) = ionic% Ye*A(i)/Z(i) * rho
                call eval_PPP_electron_table(rho_species(i),T,Z(i),K_e_species(i),ierr)
                call interp_okay% assert(ierr == 0)
            end do
            K_e = ionic% Ye/sum(Z*Y/K_e_species)
        case default
            stop 'bad selection'
        end select
    end function get_tabulated_conductivity
    subroutine write_table(prefix,tab)
        character(len=*), intent(in) :: prefix
        real(dp), dimension(NlgT,Nlgrho), intent(in) :: tab
        integer :: ierr
        character(len=16) :: method_id
        character(len=64) :: filename, head_fmt, body_fmt
        integer :: iounit, i
        type(assertion) :: io_okay=assertion(scope='write_table')
        
        call space_to_underbar(method,method_id)
        write(filename,'(a,"_",a)') trim(prefix),trim(method_id)
        select case(prefix)
        case ('Gamma','eta_e')
            write(body_fmt,'("(f8.2,",i0,"f8.1)")') NlgT
        case('lgK_fml','lgK_tab')
            write(body_fmt,'("(",i0,"f8.3)")') NlgT+1
        end select
        open(newunit=iounit,file=filename,action='write',iostat=ierr)
        call io_okay% assert(ierr == 0)
        write(head_fmt,'("(i8,",i0,"f8.2)")') NlgT
        write(iounit,trim(head_fmt)) nuclib% Z(chem_ids(1)), lgT
        do i=1, Nlgrho
            write(iounit,trim(body_fmt)) lgrho(i), tab(:,i)
        end do
        close(iounit)
    end subroutine write_table
    function linspace(jmax,xmin,xmax)
        integer, intent(in) :: jmax
        real(dp), intent(in) :: xmin, xmax
        real(dp), dimension(jmax) :: linspace
        linspace = [ (xmin + (xmax-xmin)*real(j-1,dp)/real(jmax-1,dp), j=1,jmax) ]
    end function linspace
    subroutine space_to_underbar(str, name)
      character (len=*), intent(in) :: str
      character (len=*), intent(out) :: name
      integer :: i, len
      len = len_trim(str)
      name = ''
      do i=1,len
         if (str(i:i) == ' ') then
            name(i:i) = '_'
         else
            name(i:i) = str(i:i)
         end if
      end do
    end subroutine space_to_underbar
end program tabulate_mixture
