module abuntime
    ! reader for abuntime files

contains
    
    subroutine read_abuntime(abuntime_filename, nz, nion, rho, isos, Yion, ierr)
        
        use iso_fortran_env, only: iostat_end, error_unit
        use const_def
        use nucchem_def, only : iso_name_length
        
        character(len=*), intent(in) :: abuntime_filename
        integer, intent(out) :: nz  ! number of zones
        integer, intent(out) :: nion    ! number of species
        real(dp), intent(out), dimension(:), allocatable :: rho  ! nz
        character(len=iso_name_length), dimension(:), allocatable :: isos ! nion
        real(dp), intent(out), dimension(:,:), allocatable :: Yion   ! (nion,nz)
        integer, intent(out) :: ierr
 
        integer, parameter :: default_chunk_size = 4096
        integer :: unitno, ios
        integer :: k
        real(dp) :: time,temp,EFe,EFn
        
        open(newunit=unitno,file=trim(abuntime_filename), &
        &   action='read',status='old', iostat=ierr)
        if (failure('opening'//trim(abuntime_filename),ierr)) return
        
        print *,'reading '//abuntime_filename
        read(unitno,'(1x,i5)') nion
        print *,'nion = ',nion
        allocate(isos(nion))
        
        allocate(rho(default_chunk_size),Yion(nion,default_chunk_size), &
            & stat=ierr)
        if (failure('allocating rho,Yion',ierr)) return

        nz = 1
        do
            ! size check
            if (nz > size(rho)) then
                call realloc_abuntime_arrays(2*size(rho),ierr)
                if (failure('reallocating abuntime arrays',ierr)) return
            end if
            read(unitno,'(1x,1e18.10,2e11.3,2X,2E18.10)',iostat=ios) &
                 & time,temp,rho(nz),EFe,EFn
            if (ios == iostat_end) then
                nz = nz-1
                exit
            else if (ios /= 0) then
                print *,'abnormal return: ierr = ', ierr
                return
            end if
            read(unitno,'(1x,5(a5,1pe10.3))') (isos(k),Yion(k,nz),k=1,nion)
            read(unitno,*)
            nz = nz+1
            if (modulo(nz,100) == 0)  &
            & write (error_unit,'(a)',advance='no') '.'
        end do
        print *,'got ',nz,' zones'
        close(unitno)
        
        call realloc_abuntime_arrays(nz,ierr)
        
    contains
        subroutine realloc_abuntime_arrays(newsize,ierr)
            integer, intent(in) :: newsize
            integer, intent(out) :: ierr
            integer :: currentsize, oldsize
            real(dp), dimension(:), allocatable :: tmp_rho
            real(dp), dimension(:,:), allocatable :: tmp_Yion
            
            oldsize = size(rho)
            allocate(tmp_rho(newsize),tmp_Yion(nion,newsize),stat=ierr)
            if (ierr /= 0) return
            currentsize = min(oldsize,newsize)
            tmp_rho(1:currentsize) = rho(1:currentsize)
            tmp_Yion(:,1:currentsize) = Yion(:,1:currentsize)
            deallocate(rho,Yion)
            allocate(rho(newsize),Yion(nion,newsize),stat=ierr)
            if (ierr /= 0) return
            rho(1:newsize) = tmp_rho(1:newsize)
            Yion(:,1:newsize) = tmp_Yion(:,1:newsize)
            deallocate(tmp_rho,tmp_Yion)
        end subroutine realloc_abuntime_arrays
    end subroutine read_abuntime

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

!     subroutine do_make_crust(lgP, Yion, Xneut, charged_ids, ncharged, ion_info)
!         use const_def
!         use nucchem_def
!         use nucchem_lib
!         real(dp), intent(in), dimension(:) :: lgP
!         real(dp), intent(out), dimension(:,:) :: Yion   ! (HZ90_number, size(lgP))
!         real(dp), intent(out), dimension(:) :: Xneut    ! (size(lgP))
!         integer, intent(out), dimension(HZ90_number) :: charged_ids
!         integer, intent(out) :: ncharged
!         type(composition_info_type), dimension(:), intent(out) :: ion_info   ! (size(lgP))
!         integer, dimension(max_nnuclib) :: network_indcs
!
!         real(dp), allocatable, dimension(:,:) :: X
!         integer :: Ntab, i, indx, n_indx, indx1, indx2
!         integer, dimension(HZ90_number) :: indcs
!
!         real(dp), dimension(number_layers) :: lg_Pt
!         real(dp) :: lgP1, lgP2, width, Xsum
!
!         Ntab = size(lgP)
!         allocate(X(HZ90_number,Ntab))
!
!         lg_Pt = log10(transition_pressures)
!
!         ! set the network pointers
!         indcs = [(get_nuclide_index(HZ90_network(i)),i=1,HZ90_number)]
!
!         ! and the reverse lookup
!         network_indcs = 0
!         network_indcs(indcs) = [(i,i=1,HZ90_number)]
!
!         ! for each layer set composiiton...we'll smooth the transitions in the next step
!         X = 0.0
!         ! first layer
!         indx = network_indcs(get_nuclide_index(ion_composition(1)))
!         n_indx = network_indcs(get_nuclide_index('n'))
!         where(lgP <= lg_Pt(1))
!             X(n_indx,:) = Xn(1)
!             X(indx,:) = 1.0-Xn(1)
!         end where
!
!         do i = 2, number_layers
!             indx = network_indcs(get_nuclide_index(ion_composition(i)))
!             where(lgP > lg_Pt(i-1) .and. lgP <= lg_Pt(i))
!                 X(n_indx,:) = Xn(i)
!                 X(indx,:) = 1.0-Xn(i)
!             end where
!         end do
!
!         indx = network_indcs(get_nuclide_index(ion_composition(number_layers+1)))
!         where (lgP > lg_Pt(number_layers))
!             X(n_indx,:) = Xn(number_layers+1)
!             X(indx,:) = 1.0-Xn(number_layers+1)
!         end where
!
!         ! now smooth the transitions
!         do i = 1, number_layers
!             lgP1 = lg_Pt(i) - transition_width
!             lgP2 = lg_Pt(i) + transition_width
!             width = 2.0*transition_width
!             indx1 = network_indcs(get_nuclide_index(ion_composition(i)))
!             indx2 = network_indcs(get_nuclide_index(ion_composition(i+1)))
!             where(lgP >= lgP1 .and. lgP <= lgP2)
!                 X(n_indx,:) = (Xn(i)-Xn(i+1))*cos(0.5*pi*(lgP-lgP1)/width) + Xn(i+1)
!                 X(indx1,:) = (1.0-X(n_indx,:))*cos(0.5*pi*(lgP-lgP1)/width)
!                 X(indx2,:) = (1.0-X(n_indx,:))*(1.0 - cos(0.5*pi*(lgP-lgP1)/width))
!             end where
!         end do
!
!         ! loop over and compute composition moments
!         do i = 1, Ntab
!             call compute_composition_moments(HZ90_number,indcs,X(:,i), &
!             &   ion_info(i),Xsum,ncharged,charged_ids,Yion(:,i),exclude_neutrons=.TRUE., &
!             &   abunds_are_mass_fractions=.TRUE.)
!         end do
!         Xneut = X(n_indx,:)
!
!         deallocate(X)
!     end subroutine do_make_crust
    
end module abuntime
