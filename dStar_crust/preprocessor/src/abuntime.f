module abuntime
    ! reader for abuntime files
    use const_def, only: dp
    
    real(dp), parameter :: gravity = 1.85052e14_dp
    real(dp), parameter :: mdot = 8.8e4_dp*0.3_dp
    real(dp), parameter :: default_min_P_increment = 0.0_dp
    
contains


    subroutine read_abuntime(abuntime_filename, nz, nion, ncharged, P, rho, T, &
    &   isos, Yion, Xneut, charged_ids, ion_info, ierr, min_P_increment)
        
        use iso_fortran_env, only: iostat_end, error_unit
        use const_def
        use nucchem_def
        use nucchem_lib
        
        character(len=*), intent(in) :: abuntime_filename
        integer, intent(out) :: nz  ! number of zones
        integer, intent(out) :: nion    ! number of species
        integer, intent(out) :: ncharged
        real(dp), intent(out), dimension(:), allocatable :: P, rho ! nz
        real(dp), intent(out) :: T ! reference temperature
        character(len=iso_name_length), dimension(:), allocatable :: isos ! nion
        real(dp), intent(out), dimension(:,:), allocatable :: Yion   ! (nion,nz)
        real(dp), dimension(:), intent(out), allocatable :: Xneut
        integer, intent(out), dimension(:), allocatable :: charged_ids
        type(composition_info_type), dimension(:), intent(out), allocatable :: ion_info
        integer, intent(out) :: ierr
        real(dp), intent(in), optional :: min_P_increment
 
        integer, dimension(max_nnuclib) :: network_indcs
        integer, dimension(:), allocatable :: indcs
        real(dp) :: Xsum
        integer :: i, n_indx
 
        integer, parameter :: default_chunk_size = 4096
        integer :: unitno, ios
        integer :: k,nrejected
        real(dp), dimension(:), allocatable :: abunds
        real(dp) :: time,temp,EFe,EFn,last_P,delta_P
        
        if (present(min_P_increment)) then
            delta_P = 1.0+min_P_increment
        else
            delta_P = 1.0+default_min_P_increment
        end if
        
        open(newunit=unitno,file=trim(abuntime_filename), &
        &   action='read',status='old', iostat=ierr)
        if (failure('opening'//trim(abuntime_filename),ierr)) return
        
        write(error_unit,*) 'reading '//abuntime_filename
        read(unitno,'(1x,i5)') nion
        write(error_unit,*) 'nion = ',nion
        allocate(isos(nion),abunds(nion),charged_ids(nion))
                
        allocate(P(default_chunk_size), &
        & rho(default_chunk_size),Yion(nion,default_chunk_size), &
        & ion_info(default_chunk_size), Xneut(default_chunk_size), stat=ierr)
        if (failure('allocating P, rho, Yion, ion_info, Xneut',ierr)) return

        nz = 1
        last_P = 0.0
        nrejected = 0
        do
            ! size check
            if (nz > size(rho)) then
                call realloc_abuntime_arrays(2*size(P),ierr)
                if (failure('reallocating abuntime arrays',ierr)) return
            end if
            read(unitno,'(1x,1e18.10,2e11.3,2X,2E18.10)',iostat=ios) &
                 & time,T,rho(nz),EFe,EFn
            if (ios == iostat_end) then
                nz = nz-1
                exit
            else if (ios /= 0) then
                write(error_unit,*) 'abnormal return: ierr = ', ierr
                return
            end if
            P(nz) = gravity*mdot*time       
            read(unitno,'(1x,5(a5,1pe10.3))') (isos(k),abunds(k),k=1,nion)
            read(unitno,*)
            
            if (nz == 1) then
                ! set the network pointers
                allocate(indcs(nion))
                indcs = [(get_nuclide_index(adjustl(isos(i))),i=1,nion)]

                do i = 1, nion
                    if (indcs(i) == -1) then
                        print *,isos(i),': index not found'
                    end if
                end do

                ! find the neutron indx
                network_indcs = 0
                network_indcs(indcs) = [(i,i=1,nion)]
                n_indx = network_indcs(get_nuclide_index('n'))
            end if
            
            call compute_composition_moments(nion,indcs,abunds(:), &
            &   ion_info(nz),Xsum,ncharged,charged_ids,Yion(:,nz), &
            &   exclude_neutrons=.TRUE., abunds_are_mass_fractions=.FALSE.)
            Xneut(nz) = abunds(n_indx)
            
            ! only increment nz if P has increased
            if (P(nz) > last_P*delta_P) then
                last_P = P(nz)
                nz = nz+1
            else
                nrejected = nrejected+1
            end if
            if (modulo(nz,1000) == 0)  &
            & write (error_unit,'(a)',advance='no') '.'
        end do
        write(error_unit,'(/,a,i0,a)') 'got ',nz,' zones'
        write(error_unit,'(a,i0,a)') 'rejected ',nrejected,' zones'
        close(unitno)
        
        call realloc_abuntime_arrays(nz,ierr)
        
        ! scale temperature to Kelvin
        T = T* 1.0e9_dp
        
    contains
        subroutine realloc_abuntime_arrays(newsize,ierr)
            integer, intent(in) :: newsize
            integer, intent(out) :: ierr
            integer :: currentsize, oldsize
            real(dp), dimension(:), allocatable :: tmp_rho, tmp_P, tmp_Xneut
            real(dp), dimension(:,:), allocatable :: tmp_Yion
            type(composition_info_type), dimension(:), allocatable :: tmp_ion_info
            
            oldsize = size(rho)
            allocate(tmp_rho(newsize),tmp_Yion(nion,newsize),tmp_P(newsize), &
            &   tmp_ion_info(newsize), tmp_Xneut(newsize), stat=ierr)
            if (ierr /= 0) return
            currentsize = min(oldsize,newsize)
            tmp_rho(1:currentsize) = rho(1:currentsize)
            tmp_P(1:currentsize) = P(1:currentsize)
            tmp_Yion(:,1:currentsize) = Yion(:,1:currentsize)
            tmp_ion_info(1:currentsize) = ion_info(1:currentsize)
            tmp_Xneut(1:currentsize) = Xneut(1:currentsize)
            deallocate(rho,P,Yion,ion_info,Xneut)
            allocate(rho(newsize),P(newsize),Yion(nion,newsize), &
            &   ion_info(newsize), Xneut(newsize), stat=ierr)
            if (ierr /= 0) return
            rho(1:newsize) = tmp_rho(1:newsize)
            P(1:newsize) = tmp_P(1:newsize)
            Yion(:,1:newsize) = tmp_Yion(:,1:newsize)
            ion_info(1:newsize) = tmp_ion_info(1:newsize)
            Xneut(1:newsize) = tmp_Xneut(1:newsize)
            deallocate(tmp_rho,tmp_P,tmp_Yion,tmp_ion_info,tmp_Xneut)
        end subroutine realloc_abuntime_arrays
    end subroutine read_abuntime
        
    subroutine read_abuntime_cache(cache_filename,nz,nion,isos,lgP,Yion,ierr)
        use nucchem_def, only: iso_name_length
        character(len=*), intent(in) :: cache_filename
        integer, intent(out) :: nz,nion
        character(len=iso_name_length), intent(out), dimension(:), allocatable :: isos
        real(dp), dimension(:), intent(out), allocatable :: lgP
        real(dp), dimension(:,:), intent(out), allocatable :: Yion
        integer, intent(out) :: ierr
        integer :: unitno
        
        open(newunit=unitno,file=trim(cache_filename), &
        &   action='read',status='old',form='unformatted', iostat=ierr)
        if (failure('opening'//trim(cache_filename),ierr)) return
        
        read(unitno) nz
        read(unitno) nion
        allocate(isos(nion),lgP(nz),Yion(nion,nz),stat=ierr)
        if (failure('allocating abuntime tables',ierr)) then
            close(unitno)
            return
        end if
        read(unitno) isos
        read(unitno) lgP
        read(unitno) Yion

        close(unitno)        
    end subroutine read_abuntime_cache

    subroutine write_abuntime_cache(cache_filename,nz,nion,ncharged,isos, &
    &   charged_ids,ion_info,P,Yion,ierr)
        use nucchem_def
        
        character(len=*), intent(in) :: cache_filename
        integer, intent(in) :: nz, nion, ncharged
        character(len=iso_name_length), intent(in), dimension(:) :: isos
        integer, intent(in), dimension(:) :: charged_ids
        type(composition_info_type), intent(in), dimension(:) :: ion_info
        real(dp), intent(in), dimension(:) :: P
        real(dp), intent(in), dimension(:,:) :: Yion
        integer, intent(out) :: ierr
        integer :: unitno
        
        ierr = 0
        write(*,*) 'cache_filename is ',cache_filename
        open(newunit=unitno, file=trim(cache_filename),action='write', &
        &   form='unformatted',iostat=ierr)
        if (failure('opening '//trim(cache_filename),ierr)) return

        write(unitno) nz
        write(unitno) nion
        write(unitno) ncharged
        write(unitno) isos
        write(unitno) charged_ids
        write(unitno) ion_info
        write(unitno) P
        write(unitno) Yion
        close(unitno)
    end subroutine write_abuntime_cache
    
    subroutine check_abuntime(eos_handle, rho, T, isos, P,  &
    &   Yion, Xneut, charged_ids, ncharged, ion_info, deltaRho)    
        use const_def
        use nucchem_def
        use nucchem_lib
        use superfluid_def
        use superfluid_lib
        use dStar_eos_lib

        integer, intent(in) :: eos_handle
        real(dp), dimension(:), intent(in) :: rho
        real(dp), intent(in) :: T
        character(len=iso_name_length), dimension(:), intent(in) :: isos
        real(dp), dimension(:), intent(in) :: P
        real(dp), dimension(:,:), intent(in) :: Yion
        real(dp), dimension(:), intent(in) :: Xneut
        integer, intent(in), dimension(:) :: charged_ids
        integer, intent(in) :: ncharged
        type(composition_info_type), dimension(:), intent(in) :: ion_info
        real(dp), intent(out), dimension(:) :: deltaRho
        !
        integer :: nz, nion, i
        real(dp) :: k,Tcs(max_number_sf_types),chi,kFn,kFp
        real(dp), dimension(:), allocatable :: lgP, lgRho,lgEps
        integer :: phase
        type(crust_eos_component), dimension(num_crust_eos_components) :: &
        & eos_components
        real(dp), dimension(num_dStar_eos_results) :: res

        nz = size(rho)
        nion = size(isos)

        allocate(lgP(nz), lgRho(nz), lgEps(nz))
        lgP = log10(P)
        
        call find_densities(eos_handle, &
        &   lgP,lgRho,lgEps,Yion,ncharged,charged_ids,ion_info,T)

        deltaRho = abs(10.0_dp**lgRho/rho - 1.0_dp)
!         write(*,'(7(f9.5,tr1))') &
!         & (lgP(i),lgRho(i),log10(rho(i)), lgRho(i)-log10(rho(i)), &
!         & Xneut(i),ion_info(i)% Z, ion_info(i)%A,i=1,nz)
        deallocate(lgP,lgRho,lgEps)
    end subroutine check_abuntime

    subroutine find_densities(eos_handle,lgP,lgRho,lgEps,Yion,ncharged,charged_ids,ionic,Tref)
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
        
        Pfac = 0.25*(threepisquare)**onethird *hbar*clight*avo**(4.0*onethird)
        Ntab = size(lgP)
        imax = 20
        epsx = 1.0d-8
        epsy = 1.0d-8
        
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
            if (ierr /= 0) then
                write (*,*) 'unable to bracket root',lgP(i), x1, x3, y1, y3
                cycle
            end if
            
            lgRho(i) = safe_root_with_initial_guess(match_density,lgRho_guess,x1,x3,y1,y3, &
            &   imax,epsx,epsy,lrpar,rpar,lipar,ipar,ierr)
            if (ierr /= 0) then
                write(*,*) 'unable to converge', lgP(i), x1, x3, y1, y3
                cycle
            end if
            lgEps(i) = rpar(ncharged+14)
        end do
    end subroutine find_densities
    
    real(dp) function match_density(lgRho, dfdlgRho, lrpar, rpar, lipar, ipar, ierr)
       ! returns with ierr = 0 if was able to evaluate f and df/dx at x
       ! if df/dx not available, it is okay to set it to 0
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
       chi = nuclear_volume_fraction(rho,ionic,default_nuclear_radius)
	   kFp = 0.0_dp
	   kFn = neutron_wavenumber(rho,ionic,chi)
	   call sf_get_results(kFp,kFn,Tcs)
       call eval_crust_eos(eos_handle,rho,T,ionic,ncharged, &
       &	charged_ids,Yion,Tcs,res,phase,chi)
       Eint = res(i_lnE)
       
       lgPwant = rpar(ncharged+12)
       lgP = res(i_lnP)/ln10

       rpar(ncharged+14) = log10(rho*(1.0+dot_product(Yion(:),nuclib%mass_excess(charged_ids))/amu_n + Eint/clight2))
       rpar(lrpar) = res(i_chiRho)
       ierr = 0
       match_density = lgP - lgPwant
       deallocate(charged_ids,Yion)
    end function match_density
    
    function failure(msg,ierr)
        use, intrinsic :: iso_fortran_env, only: error_unit
        character(len=*), intent(in) :: msg
        integer, intent(in) :: ierr
        logical :: failure
        
        failure = (ierr /= 0)
        if (failure) then
            write(error_unit,*) trim(msg),': ierr = ',ierr
        end if
    end function failure
    
end module abuntime
