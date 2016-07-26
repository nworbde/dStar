program print_abuntime
    
    use iso_fortran_env, only: error_unit, output_unit
    use constants_def, only: dp 
    use nucchem_def
    use abuntime
    
    real(dp), parameter :: abundance_threshold = 0.001_dp
    real(dp), parameter :: default_min_lgP_increment = 0.0001_dp
    character(len=*), parameter ::   &
        &   cache_filename = '../data/cache/abuntime_lx2_5_cache.bin'
    character(len=128) :: line_format
    integer, dimension(:), allocatable :: above_thresh
    integer :: nz, nion, ncharged
    character(len=iso_name_length), dimension(:), allocatable :: isos
    integer, dimension(:), allocatable :: charged_ids
    type(composition_info_type), dimension(:), allocatable :: ion_info
    real(dp) :: T, lgP_last
    real(dp), dimension(:), allocatable :: Xn, lgP, lgRho, lgEps
    real(dp), dimension(:,:), allocatable :: Yion
    integer :: ierr, i, k
    character(len=iso_name_length), dimension(:), allocatable :: abund_isos
    real(dp), dimension(:), allocatable :: abund_Ys
    integer :: nabund
    logical, dimension(:,:), allocatable :: threshold
    integer :: summary_unit, iso_unit
    
    call read_abuntime_cache(cache_filename,nz,nion,ncharged,isos, &
        &   charged_ids,ion_info,Xn,T,lgP,lgRho,lgEps,Yion,ierr)
    
    allocate(threshold(nion,nz),above_thresh(nz))
    threshold = Yion > abundance_threshold*maxval(Yion(:,:))
    above_thresh = count(threshold,dim=1)
    nabund = maxval(above_thresh)
    write(error_unit,'(/,a,es11.4,a,i5)')  &
    &   'max number of isotopes with Y > ', &
    &   abundance_threshold,'*Ymax is',nabund
    
    allocate(abund_isos(nabund),abund_Ys(nabund))
    write(line_format,'(a,i0,a)') "(4(f8.4,tr1),tr2,",nabund,"(a6,tr1,es11.4))"
    lgP_last = -100.0
    
    open(newunit=iso_unit,file='abuntime_isos',action='write',status='unknown')
    write(iso_unit,'(a5,tr1)') isos
    
    open(newunit=summary_unit,file='abuntime_summary',action='write', &
    &   status='unknown')
    do i = 1, nz
        if (lgP(i) > lgP_last+default_min_lgP_increment) then
            lgP_last = lgP(i)
            abund_isos(1:above_thresh(i)) = pack(isos,threshold(:,i))
            abund_Ys(1:above_thresh(i)) = pack(Yion(:,i),threshold(:,i))
            write(summary_unit,line_format) lgP(i), lgRho(i), lgEps(i), Xn(i), &
                 &  (abund_isos(k),abund_Ys(k),k=1,above_thresh(i))
        end if
    end do
    close(summary_unit)
    deallocate(threshold,above_thresh,abund_isos,abund_Ys)
    
end program print_abuntime
