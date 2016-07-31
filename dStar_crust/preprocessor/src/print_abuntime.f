program print_abuntime
    
    use iso_fortran_env, only: error_unit, output_unit
    use constants_def
    use constants_lib
    use nucchem_def
    use nucchem_lib
    use abuntime
    
    character(len=*), parameter ::   &
        &   cache_filename = '../data/cache/abuntime_lx2_5_cache.bin'
    character(len=128) :: line_format
    integer, dimension(:), allocatable :: above_thresh
    integer :: nz, nion
    character(len=iso_name_length), dimension(:), allocatable :: isos
    real(dp) :: T, lgP_last
    real(dp), dimension(:), allocatable :: lgP, Xn
    real(dp), dimension(:,:), allocatable :: Yion
    integer :: ierr, i, k, ipos, nz_kept
    character(len=iso_name_length), dimension(:), allocatable :: abund_isos
    real(dp), dimension(:), allocatable :: abund_Ys
    integer :: nabund
    logical, dimension(:,:), allocatable :: threshold
    logical, dimension(:), allocatable :: ion_mask
    integer :: summary_unit, iso_unit, ncharged, n_indx
    
    call constants_init('',ierr)
    call nucchem_init('../../data/',ierr)
        
    call read_abuntime_cache(cache_filename,nz,nion,isos,T,lgP,Yion,ierr)

    open(newunit=iso_unit,file='abuntime_isos',action='write',status='unknown')
    write(iso_unit,'(a5)') isos
    close(iso_unit)
    
    allocate(Xn(nz))
    allocate(ion_mask(nion))
    ion_mask = .TRUE.
    do i = 1, nion
        if (adjustl(isos(i)) == 'n') then
            ion_mask(i) = .FALSE.
            n_indx = i
        end if
    end do
    Xn = Yion(n_indx,:)
    ncharged = count(ion_mask)
    do i = 1, nz
        Yion(1:ncharged,i) = pack(Yion(:,i),ion_mask)
    end do
    isos(1:ncharged) = pack(isos,ion_mask)
    
    allocate(threshold(nion,nz),above_thresh(nz))
    do i = 1, nz
        threshold(1:ncharged,i) = Yion(1:ncharged,i) > abundance_threshold*maxval(Yion(1:ncharged,i))
    end do
    above_thresh = count(threshold,dim=1)
    nabund = maxval(above_thresh)
    write(error_unit,'(/,a,es11.4,a,i5)')  &
    &   'max number of isotopes with Y > ', &
    &   abundance_threshold,'*Ymax is',nabund
    
    write(error_unit,'(a,f7.3)') 'at lgP = ',lgP(maxloc(above_thresh))
    
    allocate(abund_isos(nabund),abund_Ys(nabund))
    write(line_format,'(a,i0,a)') "(f7.3,tr2,",nabund+1,"(a6,tr1,es10.3))"
    
    open(newunit=summary_unit,file='abuntime_summary',action='write', &
    &   status='unknown')
    
    do i = 1, nz
        abund_isos(1:above_thresh(i)) = pack(isos, threshold(1:ncharged,i))
        abund_Ys(1:above_thresh(i)) = pack(Yion(1:ncharged,i),threshold(1:ncharged,i))
        write(summary_unit,line_format) lgP(i), 'n', Xn(i), &
        &   (abund_isos(k),abund_Ys(k),k=1,above_thresh(i))
    end do
    
    close(summary_unit)
    deallocate(threshold,above_thresh,abund_isos,abund_Ys)
    call nucchem_shutdown
    
end program print_abuntime
