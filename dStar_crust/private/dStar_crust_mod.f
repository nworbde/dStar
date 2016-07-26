module dStar_crust_mod
    use dStar_crust_def
    integer, parameter :: crust_filename_length=128
contains
    
    subroutine do_load_crust_table(prefix,eos_handle,Tref,ierr)
        use, intrinsic :: iso_fortran_env, only: error_unit
        character(len=*), intent(in) :: prefix
        integer, intent(in) :: eos_handle
        real(dp), intent(in) :: Tref
        integer, intent(out) :: ierr
        logical, parameter :: dbg = .FALSE.
		type(crust_table_type), pointer :: tab
        real(dp), pointer, dimension(:,:) :: lgRho_val    
        character(len=crust_filename_length) :: table_name, cache_filename
        logical :: have_cache
        integer :: unitno
        
        tab => crust_table
		! if the table is already allocated, issue a warning and scrub the table
		if (tab% is_loaded) then
			write(error_unit,'(a)') 'do_load_crust_table: overwriting already loaded table'
			call do_free_crust_table(tab)
		end if

        call generate_crust_filename(prefix,Tref,table_name)
        cache_filename = trim(crust_datadir)//'/cache/'//trim(table_name)//'.bin'
        inquire(file=cache_filename,exist=have_cache)
        if (have_cache) then
            call do_read_crust_cache(cache_filename,tab,ierr)
            if (ierr == 0) return
        end if
        
        ! if we don't have the table, or could not load it, then generatue a 
        ! new one and write to cache
        tab% nv = crust_default_number_table_points
        tab% lgP_min = crust_default_lgPmin
        tab% lgP_max = crust_default_lgPmax
        call do_generate_crust_table(prefix,eos_handle,Tref,tab)
        tab% is_loaded = .TRUE.
        
        lgRho_val(1:4,1:tab% nv) => tab% lgRho(1:4*tab% nv)
        
        ! write informative message about range of table
        if (dbg) then
            write(error_unit,'(a,2f8.3)') 'do_load_crust_table: lgNb min, max = ', &
            &   10.0**(minval(lgRho_val(1,:)-log10(amu))-39.0), 10.0**(maxval(lgRho_val(1,:)-log10(amu))-39.0)
            write(error_unit,'(t21,a,2f8.3)') 'lgP min, max = ', &
            &   tab% lgP_min, tab% lgP_max
        end if
        
        if (.not.have_cache) then
            call do_write_crust_cache(cache_filename,tab,ierr)
        end if
        
    end subroutine do_load_crust_table
    
    subroutine do_generate_crust_table(prefix,eos_handle,Tref,tab)
        use constants_def, only: avogadro
        use nucchem_def, only: composition_info_type
        use hz90
        use interp_1d_def
        use interp_1d_lib
        character(len=*), intent(in) :: prefix
        integer, intent(in) :: eos_handle
        real(dp), intent(in) :: Tref
        type(crust_table_type), pointer :: tab
		real(dp), dimension(:), pointer :: work=>null()
        real(dp), pointer, dimension(:,:) :: lgRho_val, lgEps_val
        real(dp) :: lgPmin, delta_lgP
        integer :: N, i, ierr
        integer :: ncharged
        integer, dimension(HZ90_number) :: charged_ids
        real(dp), dimension(:), allocatable :: lgP, Xneut, lgRho, lgEps
        real(dp), dimension(:,:), allocatable :: Yion
        type(composition_info_type), dimension(:), allocatable :: ion_info
            
        lgPmin = tab% lgP_min
        delta_lgP = tab% lgP_max - lgPmin
        N = tab% nv
        allocate(lgP(N), Xneut(N), lgRho(N), lgEps(N), Yion(HZ90_number,N), ion_info(N))
        
        lgP = [ (lgPmin + real(i-1,dp)*(delta_lgP)/real(N-1,dp), i = 1, N)]
        
        call do_make_crust(lgP,Yion,Xneut,charged_ids,ncharged,ion_info)        
        call find_densities(eos_handle,lgP,lgRho,lgEps, &
        	& Yion,ncharged,charged_ids,ion_info,Tref)
        
        call do_allocate_crust_table(tab, N, ierr)
        lgRho_val(1:4,1:N) => tab% lgRho(1:4*N)
        lgEps_val(1:4,1:N) => tab% lgEps(1:4*N)
        
        tab% lgP = lgP
        lgRho_val(1,:) = lgRho
        lgEps_val(1,:) = lgEps
                
        allocate(work(tab% nv*pm_work_size))
        call interp_pm(tab% lgP, tab% nv, tab% lgRho, pm_work_size, work, &
        &   'do_generate_crust_table: Rho', ierr)
        call interp_pm(tab% lgP, tab% nv, tab% lgEps, pm_work_size, work, &
        &   'do_generate_crust_table: Eps', ierr)
        deallocate(work)
        deallocate(lgP, Xneut, lgRho, lgEps, Yion, ion_info)

    end subroutine do_generate_crust_table
    
    subroutine do_read_crust_cache(cache_filename,tab,ierr)
        character(len=*), intent(in) :: cache_filename
        type(crust_table_type), pointer :: tab
        integer, intent(out) :: ierr
        integer :: unitno, n
        
        open(newunit=unitno,file=trim(cache_filename), &
        &   action='read',status='old',form='unformatted', iostat=ierr)
        if (failure('opening'//trim(cache_filename),ierr)) return
        
        read(unitno) n
        call do_allocate_crust_table(tab,n,ierr)
        if (failure('allocating table',ierr)) then
            close(unitno)
            return
        end if
        
        read(unitno) tab% lgP(:)
        read(unitno) tab% lgRho
        read(unitno) tab% lgEps
        close(unitno)
        
        tab% lgP_min = minval(tab% lgP)
        tab% lgP_max = maxval(tab% lgP)
        tab% is_loaded = .TRUE.
    end subroutine do_read_crust_cache
    
    subroutine do_write_crust_cache(cache_filename,tab,ierr)
        character(len=*), intent(in) :: cache_filename
        type(crust_table_type), pointer :: tab
        integer, intent(out) :: ierr
        integer :: unitno
        
        ierr = 0
        if (tab% nv == 0) return
        open(newunit=unitno, file=trim(cache_filename),action='write', &
        &   form='unformatted',iostat=ierr)
        if (failure('opening '//trim(cache_filename),ierr)) return
        
        write(unitno) tab% nv
        write(unitno) tab% lgP
        write(unitno) tab% lgRho
        write(unitno) tab% lgEps
        close(unitno)
    end subroutine do_write_crust_cache
    
    subroutine do_allocate_crust_table(tab,n,ierr)
        type(crust_table_type),pointer :: tab
        integer, intent(in) :: n
        integer, intent(out) :: ierr
                
		allocate(tab% lgP(n), tab% lgRho(4*n), tab% lgEps(4*n), stat=ierr)
        if (ierr /= 0) return
        tab% nv = n
    end subroutine do_allocate_crust_table
    
	subroutine generate_crust_filename(prefix,Tref,filename)
		! naming convention for flies is prefix_ttt
		! where ttt = 100*log10(Tref), to 3 significant digits
		character(len=*), intent(in) :: prefix
		real(dp), intent(in) :: Tref
		character(len=crust_filename_length), intent(out) :: filename
		
		write (filename,'(a,"_",i0.3)') trim(prefix), int(100.0*log10(Tref))
	end subroutine generate_crust_filename

	subroutine do_free_crust_table(tab)
		type(crust_table_type), pointer :: tab
		tab% nv = 0
		tab% lgP_min = 0.0
		tab% lgP_max = 0.0
		if (allocated(tab% lgP)) deallocate(tab% lgP)
		if (associated(tab% lgRho)) deallocate(tab% lgRho)
		if (associated(tab% lgEps)) deallocate(tab% lgEps)
        nullify(tab% lgRho)
        nullify(tab% lgEps)
		tab% is_loaded = .FALSE.
	end subroutine do_free_crust_table

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
end module dStar_crust_mod
