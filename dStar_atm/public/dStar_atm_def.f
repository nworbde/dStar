module dStar_atm_def
    use constants_def
    
    implicit none
    integer, parameter :: atm_NT = 256
    integer, parameter :: atm_Ng = 8
    real(dp), parameter :: atm_lggrav_min = 13.5
    real(dp), parameter :: atm_lggrav_max = 15.0
    real(dp), parameter :: atm_lgTb_min = 6.5
    real(dp), parameter :: atm_lgTb_max = 9.5
    real(dp), parameter :: atm_lgTeff_min = 5.3
    real(dp), parameter :: atm_lgTeff_max = 7.0

    type atm_table_type
        logical :: is_loaded = .FALSE.
        integer :: linear_lgT
        integer :: linear_lggrav
        real(dp) :: Tb_min
        real(dp) :: Tb_max
        real(dp) :: grav_max
        real(dp) :: grav_max
        real(dp), dimension(atm_NT) :: lgTb
        real(dp), dimension(atm_Ng) :: lggrav
        real(dp), dimension(atm_NT,atm_Ng) :: lgTeff
        real(dp), dimension(atm_NT,atm_Ng) :: lgflux
    end type atm_table_type

    type(atm_table_type), target, save :: atm_table
    logical, save :: atm_is_initialized = .FALSE.
    character(len=256), save :: atm_datadir
    
end module dStar_atm_def
