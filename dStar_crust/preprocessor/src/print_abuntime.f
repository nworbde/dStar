program print_abuntime
    use exceptions_lib
    use const_def, only: dp
    use nucchem_def
    use nucchem_lib
    use composition_handler
    
    character(len=*), parameter :: datadir = '../data/'
    integer :: argument_count
    character(len=64) :: argument
    real(dp) :: abundance_threshold
    character(len=64) abuntime_stem, cache_filename, iso_filename, summary_filename
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

    character(len=128) :: alert_msg
    type(assertion) :: command_arguments=assertion(scope='print_abuntime', &
        & message='USAGE: print_composition <file stem> <abundance threshold>')
    type(assertion) :: nucchem_load_okay=assertion(scope='process_abuntime', &
        & message='unable to initialize nucchem')
    type(assertion) :: read_cache_okay=assertion(scope='print_abuntime', &
        & message='read_abuntime_cache')
    type(alert) :: status = alert(scope='print_abuntime')

    argument_count = command_argument_count()
    call command_arguments% assert(argument_count==2)
    
    call get_command_argument(1,argument)
    read(argument,*) abuntime_stem
    call get_command_argument(2,argument)
    read(argument,*) abundance_threshold
    
    cache_filename = datadir//trim(abuntime_stem)//'.bin'
    iso_filename = datadir//trim(abuntime_stem)//'_isos'
    summary_filename = datadir//trim(abuntime_stem)//'_summary'
    
    call nucchem_init('../../data/',ierr)
    call nucchem_load_okay% assert(ierr==0)
    
    call read_composition_cache(cache_filename,nz,nion,isos,T,lgP,Yion,ierr)
    call read_cache_okay% assert(ierr==0)
    open(newunit=iso_unit,file=iso_filename,action='write',status='unknown')
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
    write(alert_msg,'(a,es11.4,a,i5,a,f7.3)')  &
    &   'max number of isotopes with Y > ', &
    &   abundance_threshold,'*Ymax is',nabund,' at lgP = ',lgP(maxloc(above_thresh))
    call status% report(alert_msg)

    allocate(abund_isos(nabund),abund_Ys(nabund))
    write(line_format,'(a,i0,a)') "(f7.3,tr2,",nabund+1,"(a6,tr1,es10.3))"

    open(newunit=summary_unit,file=summary_filename,action='write', &
    &   status='unknown')

    do i = 1, nz
        abund_isos(1:above_thresh(i)) = pack(isos, threshold(1:ncharged,i))
        abund_Ys(1:above_thresh(i)) = pack(Yion(1:ncharged,i),threshold(1:ncharged,i))
        write(summary_unit,line_format) lgP(i), 'n', Xn(i), &
        &   (abund_isos(k),abund_Ys(k),k=1,above_thresh(i))
    end do

    close(summary_unit)
    deallocate(threshold,above_thresh,abund_isos,abund_Ys)
    deallocate(isos,Yion,lgP)
    call nucchem_shutdown
    
end program print_abuntime
