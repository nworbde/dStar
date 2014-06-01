module dStar_eos_def
    use constants_def, only: dp
    
    ! it can be interesting to look at the individual components; supply a vector of this to get the individual pieces
    type crust_eos_component
        real(dp) :: P
        real(dp) :: E
        real(dp) :: S
        real(dp) :: F
        real(dp) :: Cv
        ! divide the following by P, if P /= 0, to get chiRho and chiT.
        real(dp) :: dPdlnRho
        real(dp) :: dPdlnT
    end type crust_eos_component

    real(dp), parameter :: default_Gamma_melt = 175.0
    real(dp), parameter :: default_rsi_melt = 140.0
    real(dp), parameter :: default_Y_threshold = 1.0e-9
    ! warning codes
    real(dp), parameter :: Q2_threshold = 18.0
    integer, parameter :: strong_quantum_effects = 1
    
    ! for specifying the phase, liquid or solid
    integer, parameter :: liquid_phase = 1
    integer, parameter :: solid_phase = 2

    ! nuclear information: if chi, the fraction of the W-S cell occupied by a nucleus, 
    ! is set to use_default_nuclear_size, then it will be computed as 
    ! chi = (4*pi/3)A(1.13 fm)**3 * n_nuclei
    real(dp), parameter :: default_nuclear_radius = 1.13 ! fm
    real(dp), parameter :: use_default_nuclear_size = -1.0

    ! for storing the helmholtz (electron-positron) EOS.. this is adapted from MESA
    
    include 'helm_results.f'
    
    type Helm_Table
        ! controls
        logical :: with_coulomb_corrections

        ! sizes of the arrays
        integer :: imax
        integer :: jmax

        !..density and temperature
        double precision :: logtlo ! log10 temp
        double precision :: logthi ! log10 temp
        double precision :: templo ! 10**logtlo
        double precision :: temphi ! 10**logthi
        double precision :: logtstp
        double precision :: logtstpi
        double precision :: logdlo ! log10 rho
        double precision :: logdhi ! log10 rho
        double precision :: denlo ! 10**logdlo
        double precision :: denhi ! 10**logdhi
        double precision :: logdstp
        double precision :: logdstpi
        double precision, dimension(:), pointer :: d ! (imax) 
        double precision, dimension(:), pointer :: t ! (jmax) 

        !..for the helmholtz free energy tables
        double precision, dimension(:,:), pointer :: f ! (imax,jmax) 
        double precision, dimension(:,:), pointer :: fd ! (imax,jmax)
        double precision, dimension(:,:), pointer :: ft ! (imax,jmax)
        double precision, dimension(:,:), pointer :: fdd ! (imax,jmax)
        double precision, dimension(:,:), pointer :: ftt ! (imax,jmax)
        double precision, dimension(:,:), pointer :: fdt ! (imax,jmax)
        double precision, dimension(:,:), pointer :: fddt ! (imax,jmax)
        double precision, dimension(:,:), pointer :: fdtt ! (imax,jmax)
        double precision, dimension(:,:), pointer :: fddtt ! (imax,jmax)

        !..for the pressure derivative with density tables
        double precision, dimension(:,:), pointer :: dpdf ! (imax,jmax)
        double precision, dimension(:,:), pointer :: dpdfd ! (imax,jmax)
        double precision, dimension(:,:), pointer :: dpdft ! (imax,jmax)
        double precision, dimension(:,:), pointer :: dpdfdt ! (imax,jmax)

        !..for chemical potential tables
        double precision, dimension(:,:), pointer :: ef ! (imax,jmax)
        double precision, dimension(:,:), pointer :: efd ! (imax,jmax)
        double precision, dimension(:,:), pointer :: eft ! (imax,jmax)
        double precision, dimension(:,:), pointer :: efdt ! (imax,jmax)

        !..for the number density tables
        double precision, dimension(:,:), pointer :: xf ! (imax,jmax)
        double precision, dimension(:,:), pointer :: xfd ! (imax,jmax)
        double precision, dimension(:,:), pointer :: xft ! (imax,jmax)
        double precision, dimension(:,:), pointer :: xfdt ! (imax,jmax)

        !..for storing the differences
        double precision, dimension(:), pointer :: dt_sav ! (jmax)
        double precision, dimension(:), pointer :: dt2_sav ! (jmax)
        double precision, dimension(:), pointer :: dti_sav ! (jmax)
        double precision, dimension(:), pointer :: dt2i_sav ! (jmax)
        double precision, dimension(:), pointer :: dt3i_sav ! (jmax)
        double precision, dimension(:), pointer :: dd_sav ! (imax)
        double precision, dimension(:), pointer :: dd2_sav ! (imax)
        double precision, dimension(:), pointer :: ddi_sav ! (imax)
        double precision, dimension(:), pointer :: dd2i_sav ! (imax)
        double precision, dimension(:), pointer :: dd3i_sav ! (jmax)

    end type Helm_Table

    type (Helm_Table), pointer :: crust_eos_ht
    type (Helm_Table), pointer :: eos_ht
    
end module dStar_eos_def
