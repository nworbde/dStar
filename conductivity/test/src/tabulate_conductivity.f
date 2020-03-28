program tabulate_conductivity
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
    ! set the table limits and resolution in lg[rho/(g cm**-3)]
    integer :: Nlgrho
    real(dp), parameter :: lgrho_min=5.0_dp, lgrho_max = 9.7_dp, delta_lgrho = 0.1_dp
    ! set the table limits and resolution in lg[T/K]
    integer :: NlgT
    real(dp), parameter :: lgT_min=6.5_dp, lgT_max = 9.0_dp, delta_lgT = 0.1_dp
    ! storage for lg(rho), lg(T)
    real(dp), dimension(:), allocatable :: lgrho, lgT
    ! storage for lg(K) (conductivity, both fml and table), eta_e (eletron degeneracy),
    ! Gamma (Coulomb pot. between ions/kT) 
    real(dp), dimension(:,:), allocatable :: lgK_fml, lgK_tab, eta_e, Gamma

    ! user inputs
    ! directory for EOS and composition datafiles
    character(len=128) :: datadir
    ! name of nuclide
    character(len=iso_name_length) :: nuclide
    
    ! for interfacing with eos, chem, and conductivity modules
    ! for the nuclide database
    integer, dimension(1) :: chem_ids, charged_ids
    real(dp), dimension(1) :: Y, Yion
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
    character(len=64) :: progname
    ! miscellaneous
    integer :: i,j,k
    real(dp) :: rho,T,K_e,eta,mu_e
    
    ! get user inputs
    if (command_argument_count() /= 2) then
        call get_command_argument(0,progname)
        stop 'Usage: '//trim(progname)//' <path/of/dStar/data/directory> <species>'
    end if
    call get_command_argument(1,datadir)
    call get_command_argument(2,nuclide)
    
    ! turn off all status messages except for warnings and failures
    call set_verbosity(1)
    
    ! initialization
    ! physical constants
    call constants_init('',ierr)
    call check_okay% assert(ierr==0)
    
    ! isotope data
    call nucchem_init(datadir,ierr)
    call check_okay% assert(ierr==0)
    
    ! equation of state
    call dStar_eos_startup(datadir)
    eos_handle = alloc_dStar_eos_handle(ierr)
    call check_okay% assert(ierr==0)
    ! suppress the warnings for degenerate ions
    call dStar_eos_set_controls(eos_handle,suppress_warnings=.TRUE.)
    
    ! thermal conductivity
    call conductivity_startup(datadir)
    cond_handle = alloc_conductivity_handle(ierr)
    call check_okay% assert(ierr==0)
    ! set electron-electron scattering formula to Potekhin, Chabrier and Yakovlev (1997)
    call conductivity_set_controls(cond_handle,which_ee_scattering=icond_pcy)

    ! set superfluid temps to a uniform value; we aren't using these, so the precise value is unimportant
    Tcs = 1.0e9_dp

    ! set the chemistry
    call check_okay% set_message('locating nuclide information')
    chem_ids = [ get_nuclide_index(StrLowCase(nuclide)) ]
    call check_okay% assert(chem_ids(1) /= nuclide_not_found)
    
    Y = 1.0_dp/nuclib% A(chem_ids(1))
    call compute_composition_moments(1,chem_ids,Y,ionic,Xsum, &
        & ncharged, charged_ids, Yion, exclude_neutrons=.TRUE.)

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
            call eval_PPP_electron_table(rho,T,ionic% Z,K_e,ierr)
            lgK_tab(j,i) = log10(K_e)
        end do temperature
    end do density
    
    ! write out tables for lg(K), mu_e, eta, Gamma
    call write_table('lgK_fml',lgK_fml)
    call write_table('eta_e',eta_e)
    call write_table('Gamma',Gamma)
    call write_table('lgK_tab',lgK_tab)

    call conductivity_shutdown
    call dStar_eos_shutdown
    call nucchem_shutdown
    
contains
    subroutine write_table(prefix,tab)
        character(len=*), intent(in) :: prefix
        real(dp), dimension(NlgT,Nlgrho), intent(in) :: tab
        integer :: ierr
        character(len=64) :: filename, head_fmt, body_fmt
        integer :: iounit, i
        type(assertion) :: io_okay=assertion(scope='write_table')
        
        write(filename,'(a,"_",a)') trim(prefix),trim(nuclib% name(chem_ids(1)))
        select case(prefix)
        case ('Gamma','eta_e')
            write(body_fmt,'("(f8.2,",i0,"f8.1)")') NlgT
        case('lgK_fml','lgK_tab')
            write(body_fmt,'("(",i0,"f8.2)")') NlgT+1
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
end program tabulate_conductivity
