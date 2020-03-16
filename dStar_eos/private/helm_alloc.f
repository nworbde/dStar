! This is modified from the MESA file
! The header for MESA is as follows.
!
c ***********************************************************************
!
!   Copyright (C) 2006, 2007  Bill Paxton, Frank Timmes
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
c ***********************************************************************

      module helm_alloc

      
      implicit none


      contains
            
      
      subroutine alloc_helm_table(h, imax, jmax, ierr)
         ! This routine allocates a Helm_Table and places pointer to it in h.
         ! It also allocates the arrays in the Helm_Table record.
         
         use, intrinsic :: iso_fortran_env, only: error_unit
         use exceptions_lib
         use dStar_eos_def
         
         type (Helm_Table), pointer :: h
         integer, intent(in) :: imax, jmax
         integer, intent(out) :: ierr ! 0 means AOK.
         type(failure) :: allocation_error=failure(scope='alloc_helm_table')

         ierr = 0
         
         allocate(h,stat=ierr)
         if (allocation_error% raised(ierr,'failure in attempt to allocate Helm_Table storage')) return
         
         h% imax = imax
         h% jmax = jmax 
         h% with_coulomb_corrections = .true.
         
         call alloc_1d_array(h% d, imax)
         call alloc_1d_array(h% t, jmax)
         
         call alloc_2d_array(h% f, imax, jmax)
         call alloc_2d_array(h% fd, imax, jmax)
         call alloc_2d_array(h% ft, imax, jmax)
         call alloc_2d_array(h% fdd, imax, jmax)
         call alloc_2d_array(h% ftt, imax, jmax)
         call alloc_2d_array(h% fdt, imax, jmax)
         call alloc_2d_array(h% fddt, imax, jmax)
         call alloc_2d_array(h% fdtt, imax, jmax)
         call alloc_2d_array(h% fddtt, imax, jmax)

         !..for the pressure derivative with density tables
         call alloc_2d_array(h% dpdf, imax, jmax)
         call alloc_2d_array(h% dpdfd, imax, jmax)
         call alloc_2d_array(h% dpdft, imax, jmax)
         call alloc_2d_array(h% dpdfdt, imax, jmax)

         !..for chemical potential tables
         call alloc_2d_array(h% ef, imax, jmax)
         call alloc_2d_array(h% efd, imax, jmax)
         call alloc_2d_array(h% eft, imax, jmax)
         call alloc_2d_array(h% efdt, imax, jmax)

         !..for the number density tables
         call alloc_2d_array(h% xf, imax, jmax)
         call alloc_2d_array(h% xfd, imax, jmax)
         call alloc_2d_array(h% xft, imax, jmax)
         call alloc_2d_array(h% xfdt, imax, jmax)

         !..for storing the differences
         call alloc_1d_array(h% dt_sav, jmax)
         call alloc_1d_array(h% dt2_sav, jmax)
         call alloc_1d_array(h% dti_sav, jmax)
         call alloc_1d_array(h% dt2i_sav, jmax)
         call alloc_1d_array(h% dt3i_sav, jmax)
         
         call alloc_1d_array(h% dd_sav, imax)
         call alloc_1d_array(h% dd2_sav, imax)
         call alloc_1d_array(h% ddi_sav, imax)
         call alloc_1d_array(h% dd2i_sav, imax)
         call alloc_1d_array(h% dd3i_sav, imax)

         contains
         
         subroutine alloc_1d_array(ptr,sz)
            double precision, dimension(:), pointer :: ptr
            integer, intent(in) :: sz        
            allocate(ptr(sz),stat=ierr)
            if (allocation_error% raised(ierr,'failure in attempt to allocate Helm_Table storage')) return
         end subroutine alloc_1d_array
         
         subroutine alloc_2d_array(ptr,sz1,sz2)
            double precision, dimension(:,:), pointer :: ptr
            integer, intent(in) :: sz1,sz2         
            allocate(ptr(sz1,sz2),stat=ierr)
            if (allocation_error% raised(ierr,'failure in attempt to allocate Helm_Table storage')) return
         end subroutine alloc_2d_array
      
      
      end subroutine alloc_helm_table


      subroutine read_helm_table(h, data_dir, ierr)
         use iso_fortran_env, only : error_unit
         use exceptions_lib
         use dStar_eos_def

         implicit none

      type (Helm_Table), pointer :: h
      character(*), intent(IN) :: data_dir
      integer, intent(out) :: ierr
      type(alert) :: status=alert(scope='read_helm_table')
      type(failure) :: io_failure=failure(scope='read_helm_table')

!..this routine reads the helmholtz eos file, and 
!..must be called once before the helmeos routine is invoked.

!..declare local variables
      character (len=256) :: filename, message
      integer          i,j,ios,imax,jmax
      double precision tsav,dsav,dth,dt2,dti,dt2i,dt3i,
     1                 dd,dd2,ddi,dd2i,dd3i
      
       ierr = 0

!..read the normal helmholtz free energy table
!       tlo   = 4.0d0
!       thi   = 11.0d0
!       dlo   = -10.0d0
!       dhi   = 11.0d0

!..for the bigger table
       h% logtlo   = 3.0d0
       h% logthi   = 13.0d0
       h% logdlo   = -12.0d0
       h% logdhi   = 14.0d0
       
       h% templo = 10**h% logtlo
       h% temphi = 10**h% logthi
       h% denlo = 10**h% logdlo
       h% denhi = 10**h% logdhi

       imax = h% imax
       jmax = h% jmax
       h% logtstp  = (h% logthi - h% logtlo)/float(jmax-1)
       h% logtstpi = 1.0d0/h% logtstp
       h% logdstp  = (h% logdhi - h% logdlo)/float(imax-1)
       h% logdstpi = 1.0d0/h% logdstp

       write(filename,'(2a)') trim(data_dir), '/cache/helm_table.bin'
      
       open(unit=19,file=trim(filename),
     >         action='read',status='old',iostat=ios,form='unformatted')
         
       if (ios .eq. 0) then
         
          read(19) imax
          read(19) jmax
         
         if (imax /= h% imax .or. jmax /= h% jmax) then
            ios = 1 ! wrong cached info
         else
             read(19) h% f(1:imax,1:jmax)
             read(19) h% fd(1:imax,1:jmax)
             read(19) h% ft(1:imax,1:jmax)
             read(19) h% fdd(1:imax,1:jmax)
             read(19) h% ftt(1:imax,1:jmax)
             read(19) h% fdt(1:imax,1:jmax)
             read(19) h% fddt(1:imax,1:jmax)
             read(19) h% fdtt(1:imax,1:jmax)
             read(19) h% fddtt(1:imax,1:jmax)
             read(19) h% dpdf(1:imax,1:jmax)
             read(19) h% dpdfd(1:imax,1:jmax)
             read(19) h% dpdft(1:imax,1:jmax)
             read(19) h% dpdfdt(1:imax,1:jmax)
             read(19) h% ef(1:imax,1:jmax)
             read(19) h% efd(1:imax,1:jmax)
             read(19) h% eft(1:imax,1:jmax)
             read(19) h% efdt(1:imax,1:jmax)
             read(19) h% xf(1:imax,1:jmax)
             read(19) h% xfd(1:imax,1:jmax)
             read(19) h% xft(1:imax,1:jmax)
             read(19) h% xfdt(1:imax,1:jmax)
         
            do j=1,jmax
               tsav = h% logtlo + (j-1)*h% logtstp
               h% t(j) = 10.0d0**(tsav)
            enddo
            do i=1,imax
               dsav = h% logdlo + (i-1)*h% logdstp
               h% d(i) = 10.0d0**(dsav)
            enddo
         end if
         
         close(unit=19)
         
       end if

       if (ios .ne. 0) then
      
          write(filename,'(2a)') trim(data_dir), '/helm_table.dat'
          call status% report('read '//trim(filename))
          
          ios = 0
          open(unit=19,file=trim(filename),action='read',status='old',iostat=ios)
          if (io_failure% raised(ios,'failed to open '//trim(filename))) then 
            ierr = -1
            return
          end if
         
          do j=1,jmax
            tsav = h% logtlo + (j-1)*h% logtstp
            h% t(j) = 10.0d0**(tsav)
            do i=1,imax
               dsav = h% logdlo + (i-1)*h% logdstp
               h% d(i) = 10.0d0**(dsav)
               read(19,*) 
     >            h% f(i,j), h% fd(i,j), h% ft(i,j),
     >            h% fdd(i,j), h% ftt(i,j), h% fdt(i,j),
     >            h% fddt(i,j), h% fdtt(i,j), h% fddtt(i,j)
            enddo
          enddo

         !..read the pressure derivative with density table
          do j=1,jmax
           do i=1,imax
            read(19,*) h% dpdf(i,j), h% dpdfd(i,j), h% dpdft(i,j), h% dpdfdt(i,j)
           enddo
          enddo

         !..read the electron chemical potential table
          do j=1,jmax
           do i=1,imax
            read(19,*) h% ef(i,j), h% efd(i,j), h% eft(i,j), h% efdt(i,j)
           enddo
          enddo

         !..read the number density table
          do j=1,jmax
           do i=1,imax
            read(19,*) h% xf(i,j), h% xfd(i,j), h% xft(i,j), h% xfdt(i,j)
           enddo
          enddo

          close(unit=19)
          !..write cachefile
      
          write(filename,'(2a)') trim(data_dir), '/cache/helm_table.bin'
          call status% report('write '//trim(filename))
          open(unit=19,file=trim(filename),status='replace',
     >            iostat=ios,action='write',form='unformatted')
         
          if (ios == 0) then
      
             write(19) imax
             write(19) jmax
             write(19) h% f(1:imax,1:jmax)
             write(19) h% fd(1:imax,1:jmax)
             write(19) h% ft(1:imax,1:jmax)
             write(19) h% fdd(1:imax,1:jmax)
             write(19) h% ftt(1:imax,1:jmax)
             write(19) h% fdt(1:imax,1:jmax)
             write(19) h% fddt(1:imax,1:jmax)
             write(19) h% fdtt(1:imax,1:jmax)
             write(19) h% fddtt(1:imax,1:jmax)
             write(19) h% dpdf(1:imax,1:jmax)
             write(19) h% dpdfd(1:imax,1:jmax)
             write(19) h% dpdft(1:imax,1:jmax)
             write(19) h% dpdfdt(1:imax,1:jmax)
             write(19) h% ef(1:imax,1:jmax)
             write(19) h% efd(1:imax,1:jmax)
             write(19) h% eft(1:imax,1:jmax)
             write(19) h% efdt(1:imax,1:jmax)
             write(19) h% xf(1:imax,1:jmax)
             write(19) h% xfd(1:imax,1:jmax)
             write(19) h% xft(1:imax,1:jmax)
             write(19) h% xfdt(1:imax,1:jmax)
             close(unit=19)

          end if
            
       
       end if 

!..construct the temperature and density deltas and their inverses 
       do j=1,jmax-1
        dth         = h% t(j+1) - h% t(j)
        dt2         = dth * dth
        dti         = 1.0d0/dth
        dt2i        = 1.0d0/dt2
        dt3i        = dt2i*dti
        h% dt_sav(j)   = dth
        h% dt2_sav(j)  = dt2
        h% dti_sav(j)  = dti
        h% dt2i_sav(j) = dt2i
        h% dt3i_sav(j) = dt3i
       end do
       do i=1,imax-1
        dd          = h% d(i+1) - h% d(i)
        dd2         = dd * dd
        ddi         = 1.0d0/dd
        dd2i        = 1.0d0/dd2
        dd3i        = dd2i*ddi
        h% dd_sav(i)   = dd
        h% dd2_sav(i)  = dd2
        h% ddi_sav(i)  = ddi
        h% dd2i_sav(i) = dd2i
        h% dd3i_sav(i) = dd3i
       enddo

      return

      end subroutine read_helm_table




      subroutine free_helm_table(h)
         use dStar_eos_def
         
         type (Helm_Table), pointer :: h
         
         call do_free(h% d)
         call do_free(h% t)
         
         call do_free2(h% f)
         call do_free2(h% fd)
         call do_free2(h% ft)
         call do_free2(h% fdd)
         call do_free2(h% ftt)
         call do_free2(h% fdt)
         call do_free2(h% fddt)
         call do_free2(h% fdtt)
         call do_free2(h% fddtt)

         !..for the pressure derivative with density tables
         call do_free2(h% dpdf)
         call do_free2(h% dpdfd)
         call do_free2(h% dpdft)
         call do_free2(h% dpdfdt)

         !..for chemical potential tables
         call do_free2(h% ef)
         call do_free2(h% efd)
         call do_free2(h% eft)
         call do_free2(h% efdt)

         !..for the number density tables
         call do_free2(h% xf)
         call do_free2(h% xfd)
         call do_free2(h% xft)
         call do_free2(h% xfdt)

         !..for storing the differences
         call do_free(h% dt_sav)
         call do_free(h% dt2_sav)
         call do_free(h% dti_sav)
         call do_free(h% dt2i_sav)
         call do_free(h% dt3i_sav)
         call do_free(h% dd_sav)
         call do_free(h% dd2_sav)
         call do_free(h% ddi_sav)
         call do_free(h% dd2i_sav)
         call do_free(h% dd3i_sav)
         
         deallocate(h)
         nullify(h)
         
         contains
         
         subroutine do_free(array_ptr)
            double precision, pointer :: array_ptr(:)
            if (associated(array_ptr)) deallocate(array_ptr)
         end subroutine do_free
         
         subroutine do_free2(array_ptr)
            double precision, pointer :: array_ptr(:,:)
            if (associated(array_ptr)) deallocate(array_ptr)
         end subroutine do_free2


      end subroutine free_helm_table
      end module helm_alloc
